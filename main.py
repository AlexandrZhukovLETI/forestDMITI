#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

try:
    from scipy.optimize import differential_evolution, minimize_scalar
except Exception:  # pragma: no cover
    differential_evolution = None
    minimize_scalar = None


ZALGALLER_LAMBDA = 2.27829164144


@dataclass
class SearchResult:
    n: int
    r: float
    length: float
    ell: float
    angles: np.ndarray
    success: bool
    message: str


def parse_float_list(text: str) -> list[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def parse_int_list(text: str) -> list[int]:
    return [int(x.strip()) for x in text.split(",") if x.strip()]


def decode_angles(z: np.ndarray, n: int, mode: str) -> np.ndarray:
    """Преобразует параметры оптимизации в направления n звеньев.

    mode='convex': ищем выпуклую монотонно поворачивающую ломаную.
                   Первое звено фиксируется под углом 0, остальные углы сортируются.
    mode='mirror': зеркально-симметричная выпуклая ломаная относительно вертикальной оси.
                   Это полезно как более узкий и быстрый класс кандидатов.
    mode='free':   свободная ломаная: z задает последовательные повороты.
    """
    z = np.asarray(z, dtype=float)

    if n == 1:
        return np.array([0.0])

    if mode == "convex":
        # Фиксация первого угла = 0 не ограничивает общность, потому что theta потом
        # перебирает глобальный поворот всей ломаной.
        return np.sort(np.r_[0.0, z[: n - 1]])

    if mode == "mirror":
        # Направления вида a_1,...,a_m, pi-a_m,...,pi-a_1.
        # Для нечетного n добавляется центральное звено под pi/2.
        half = n // 2
        if half == 0:
            return np.array([0.0])
        a = np.sort(np.r_[0.0, z[: max(0, half - 1)]])
        a = np.clip(a, 0.0, math.pi / 2)
        if n % 2 == 0:
            angles = np.r_[a, math.pi - a[::-1]]
        else:
            angles = np.r_[a, math.pi / 2, math.pi - a[::-1]]
        return angles[:n]

    if mode == "free":
        turns = z[: n - 1]
        angles = np.r_[0.0, np.cumsum(turns)]
        # Нормировка только для численной устойчивости, геометрию не меняет.
        return (angles + math.pi) % (2 * math.pi) - math.pi

    raise ValueError(f"Неизвестный режим: {mode}")


def bounds_for_mode(n: int, mode: str) -> list[tuple[float, float]]:
    if n <= 1:
        return []
    if mode == "convex":
        return [(0.0, math.pi)] * (n - 1)
    if mode == "mirror":
        return [(0.0, math.pi / 2)] * max(0, n // 2 - 1)
    if mode == "free":
        return [(-math.pi, math.pi)] * (n - 1)
    raise ValueError(f"Неизвестный режим: {mode}")


def polyline_vertices(angles: np.ndarray, ell: float = 1.0) -> np.ndarray:
    """Возвращает вершины p_0,...,p_n для n звеньев длины ell."""
    steps = ell * np.column_stack((np.cos(angles), np.sin(angles)))
    pts = np.vstack((np.zeros(2), np.cumsum(steps, axis=0)))
    return pts


def spans_for_thetas(points: np.ndarray, thetas: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Проекции повернутой ломаной на оси x,y для всех theta.

    Возвращает span_x(theta), span_y(theta): ширины bounding box.
    """
    x = points[:, 0]
    y = points[:, 1]
    c = np.cos(thetas)[:, None]
    s = np.sin(thetas)[:, None]

    qx = c * x[None, :] - s * y[None, :]
    qy = s * x[None, :] + c * y[None, :]

    span_x = qx.max(axis=1) - qx.min(axis=1)
    span_y = qy.max(axis=1) - qy.min(axis=1)
    return span_x, span_y


def need_at_theta(points_unit: np.ndarray, r: float, theta: float) -> tuple[float, float, float]:
    """Значение порога ell для одного угла theta."""
    c = math.cos(theta)
    st = math.sin(theta)
    x = points_unit[:, 0]
    y = points_unit[:, 1]
    qx = c * x - st * y
    qy = st * x + c * y
    sx = float(qx.max() - qx.min())
    sy = float(qy.max() - qy.min())
    eps = 1e-12
    tx = r / sx if sx > eps else float("inf")
    ty = 1.0 / sy if sy > eps else float("inf")
    return min(tx, ty), sx, sy


def required_ell(
    points_unit: np.ndarray,
    r: float,
    thetas: np.ndarray,
    refine_theta: bool = False,
    top_k: int = 3,
) -> tuple[float, int, float, float]:
    """Минимальная длина звена ell для данной формы ломаной.

    Для каждого theta плохой старт существует, если:
        ell*span_x(theta) <= r  и  ell*span_y(theta) <= 1.
    Значит, чтобы плохого старта не было, нужно:
        ell > min(r/span_x(theta), 1/span_y(theta))
    для каждого theta. Поэтому ell берется как максимум по theta от этого минимума.

    refine_theta=True дополнительно уточняет максимум по theta локальным одномерным поиском.
    Это защищает оптимизатор от "обмана сетки", когда худший угол попадает между узлами.
    """
    sx, sy = spans_for_thetas(points_unit, thetas)
    eps = 1e-12
    tx = np.full_like(sx, np.inf, dtype=float)
    ty = np.full_like(sy, np.inf, dtype=float)
    np.divide(r, sx, out=tx, where=sx > eps)
    np.divide(1.0, sy, out=ty, where=sy > eps)
    need_by_theta = np.minimum(tx, ty)
    idx = int(np.argmax(need_by_theta))
    best = float(need_by_theta[idx])
    best_sx = float(sx[idx])
    best_sy = float(sy[idx])
    best_idx = idx

    if refine_theta and minimize_scalar is not None and len(thetas) >= 8:
        step = float(2 * math.pi / len(thetas))
        # Уточняем несколько лучших узлов, а не только один: максимум может быть плоским
        # или рядом с точкой смены опорных вершин.
        candidates = np.argpartition(need_by_theta, -min(top_k, len(thetas)))[-min(top_k, len(thetas)):]
        for i in candidates:
            center = float(thetas[int(i)])

            def neg_need(t: float) -> float:
                val, _, _ = need_at_theta(points_unit, r, t % (2 * math.pi))
                return -val

            opt = minimize_scalar(neg_need, bounds=(center - step, center + step), method="bounded")
            val, local_sx, local_sy = need_at_theta(points_unit, r, float(opt.x) % (2 * math.pi))
            if val > best:
                best = float(val)
                best_sx = local_sx
                best_sy = local_sy
                best_idx = int(i)

    return best, best_idx, best_sx, best_sy

def objective_factory(n: int, r: float, mode: str, thetas: np.ndarray, refine_theta: bool = False):
    def objective(z: np.ndarray) -> float:
        angles = decode_angles(z, n, mode)
        pts = polyline_vertices(angles, ell=1.0)
        ell, _, _, _ = required_ell(pts, r, thetas, refine_theta=refine_theta)
        return n * ell

    return objective


def random_search(
    n: int,
    r: float,
    mode: str,
    thetas: np.ndarray,
    samples: int,
    seed: int,
) -> tuple[np.ndarray, float, str]:
    rng = np.random.default_rng(seed)
    bounds = np.array(bounds_for_mode(n, mode), dtype=float)
    dim = len(bounds)
    obj = objective_factory(n, r, mode, thetas)

    if dim == 0:
        z = np.array([])
        return z, obj(z), "n=1: параметров углов нет"

    best_z = None
    best_val = float("inf")
    lo = bounds[:, 0]
    hi = bounds[:, 1]

    for _ in range(samples):
        z = lo + rng.random(dim) * (hi - lo)
        val = obj(z)
        if val < best_val:
            best_val = val
            best_z = z.copy()

    return best_z, best_val, f"random search: {samples} случайных кандидатов"


def optimize_one(
    n: int,
    r: float,
    mode: str,
    theta_grid: int,
    generations: int,
    popsize: int,
    seed: int,
    random_samples: int,
    polish: bool = True,
    refine_theta: bool = False,
) -> SearchResult:
    thetas = np.linspace(0.0, 2 * math.pi, theta_grid, endpoint=False)
    bounds = bounds_for_mode(n, mode)
    obj = objective_factory(n, r, mode, thetas, refine_theta=refine_theta)

    if len(bounds) == 0:
        z = np.array([])
        length = obj(z)
        angles = decode_angles(z, n, mode)
        pts = polyline_vertices(angles, ell=1.0)
        ell, _, _, _ = required_ell(pts, r, thetas, refine_theta=refine_theta)
        return SearchResult(n, r, length, ell, angles, True, "n=1")

    if differential_evolution is not None:
        res = differential_evolution(
            obj,
            bounds=bounds,
            maxiter=generations,
            popsize=popsize,
            tol=1e-7,
            polish=polish,
            seed=seed,
            updating="immediate",
            workers=1,
        )
        z_best = np.asarray(res.x, dtype=float)
        msg = f"differential_evolution: {res.message}"
        success = bool(res.success)
    else:
        z_best, _, msg = random_search(n, r, mode, thetas, random_samples, seed)
        success = True

    angles = decode_angles(z_best, n, mode)
    pts = polyline_vertices(angles, ell=1.0)
    ell, _, _, _ = required_ell(pts, r, thetas, refine_theta=refine_theta)
    length = n * ell
    return SearchResult(n, r, length, ell, angles, success, msg)


def verify_on_finer_grid(result: SearchResult, theta_grid: int, refine_theta: bool = True) -> tuple[float, float, int]:
    thetas = np.linspace(0.0, 2 * math.pi, theta_grid, endpoint=False)
    pts_unit = polyline_vertices(result.angles, ell=1.0)
    ell, idx, _, _ = required_ell(pts_unit, result.r, thetas, refine_theta=refine_theta)
    return result.n * ell, float(thetas[idx]), idx


def save_plot(result: SearchResult, out_png: Path) -> None:
    import matplotlib.pyplot as plt

    pts = polyline_vertices(result.angles, ell=result.ell)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(pts[:, 0], pts[:, 1], marker="o")
    for k, (x, y) in enumerate(pts):
        ax.text(x, y, f" P{k}", fontsize=9)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)
    ax.set_title(f"n={result.n}, r={result.r:g}, L≈{result.length:.6f}, ell≈{result.ell:.6f}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def save_scan_plot(csv_path: Path, out_png: Path) -> None:
    import matplotlib.pyplot as plt

    rows = []
    with csv_path.open("r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append(row)

    fig, ax = plt.subplots(figsize=(7, 5))
    n_values = sorted({int(row["n"]) for row in rows})
    for n in n_values:
        sub = [row for row in rows if int(row["n"]) == n]
        sub.sort(key=lambda row: float(row["r"]))
        xs = [float(row["r"]) for row in sub]
        ys = [float(row["length"]) for row in sub]
        ax.plot(xs, ys, marker="o", label=f"n={n}")
    ax.axhline(ZALGALLER_LAMBDA, linestyle="--", label="Залгаллер 2.2782916")
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("r = длина / ширина прямоугольника")
    ax.set_ylabel("найденная длина L")
    ax.set_title("Численный поиск в серой зоне")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def print_result(result: SearchResult, fine_theta_grid: int, plot_path: Optional[Path]) -> None:
    fine_length, worst_theta, _ = verify_on_finer_grid(result, fine_theta_grid)
    excess = (fine_length / ZALGALLER_LAMBDA - 1.0) * 100.0

    print("\n=== Результат ===")
    print(f"n = {result.n}")
    print(f"r = {result.r}")
    print(f"режим поиска: см. аргумент --mode")
    print(f"L на рабочей theta-сетке  ≈ {result.length:.9f}")
    print(f"L на тонкой theta-сетке   ≈ {fine_length:.9f}")
    print(f"ell = L/n                ≈ {fine_length / result.n:.9f}")
    print(f"сравнение с 2.27829164144: {excess:+.4f}%")
    print(f"худший угол theta ≈ {worst_theta:.6f} рад = {math.degrees(worst_theta):.3f}°")
    print("углы звеньев alpha_i в градусах:")
    print("  " + ", ".join(f"{math.degrees(a):.3f}" for a in result.angles))
    print("сообщение оптимизатора:")
    print("  " + result.message)
    if plot_path is not None:
        print(f"график ломаной сохранен в: {plot_path}")


def run_scan(args: argparse.Namespace) -> None:
    n_values = parse_int_list(args.n_values)
    r_values = parse_float_list(args.r_values)
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    total = len(n_values) * len(r_values)
    done = 0

    for n in n_values:
        for r in r_values:
            done += 1
            print(f"[{done}/{total}] поиск: n={n}, r={r}", flush=True)
            res = optimize_one(
                n=n,
                r=r,
                mode=args.mode,
                theta_grid=args.theta_grid,
                generations=args.generations,
                popsize=args.popsize,
                seed=args.seed + 1000 * n + int(100 * r),
                random_samples=args.random_samples,
                polish=not args.no_polish,
                refine_theta=args.refine_theta,
            )
            fine_length, worst_theta, _ = verify_on_finer_grid(res, args.fine_theta_grid)
            rows.append(
                {
                    "n": n,
                    "r": r,
                    "length": fine_length,
                    "ell": fine_length / n,
                    "excess_vs_zalgaller_percent": (fine_length / ZALGALLER_LAMBDA - 1.0) * 100.0,
                    "worst_theta_rad": worst_theta,
                    "angles_deg": " ".join(f"{math.degrees(a):.6f}" for a in res.angles),
                }
            )

    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nCSV сохранен: {out_csv}")
    if args.plot:
        out_png = out_csv.with_suffix(".png")
        save_scan_plot(out_csv, out_png)
        print(f"график сохранен: {out_png}")


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Поиск ломаной для выхода из прямоугольного леса")
    parser.add_argument("--n", type=int, default=6, help="число равных звеньев")
    parser.add_argument("--r", type=float, default=1.5, help="отношение сторон прямоугольника: [0,r]x[0,1]")
    parser.add_argument("--mode", choices=["convex", "mirror", "free"], default="convex", help="класс ломаных")
    parser.add_argument("--theta-grid", type=int, default=720, help="число углов theta при оптимизации")
    parser.add_argument("--fine-theta-grid", type=int, default=10000, help="тонкая проверка найденного результата")
    parser.add_argument("--generations", type=int, default=100, help="число поколений differential evolution")
    parser.add_argument("--popsize", type=int, default=12, help="множитель размера популяции")
    parser.add_argument("--seed", type=int, default=1, help="seed генератора")
    parser.add_argument("--random-samples", type=int, default=20000, help="fallback, если scipy недоступен")
    parser.add_argument("--no-polish", action="store_true", help="отключить финальную локальную доводку scipy")
    parser.add_argument("--refine-theta", action="store_true", help="уточнять худший угол theta прямо внутри оптимизации; медленнее, но надежнее")
    parser.add_argument("--plot", action="store_true", help="сохранить график")
    parser.add_argument("--out-png", default="forest_polyline.png", help="имя PNG для одиночного запуска")

    parser.add_argument("--scan", action="store_true", help="запустить серию по нескольким n и r")
    parser.add_argument("--n-values", default="2,4,6,8", help="n для scan, через запятую")
    parser.add_argument("--r-values", default="1,1.25,1.5,2,3,5", help="r для scan, через запятую")
    parser.add_argument("--out-csv", default="forest_scan.csv", help="CSV для scan")

    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.scan:
        run_scan(args)
        return 0

    result = optimize_one(
        n=args.n,
        r=args.r,
        mode=args.mode,
        theta_grid=args.theta_grid,
        generations=args.generations,
        popsize=args.popsize,
        seed=args.seed,
        random_samples=args.random_samples,
        polish=not args.no_polish,
        refine_theta=args.refine_theta,
    )

    plot_path = None
    if args.plot:
        plot_path = Path(args.out_png)
        save_plot(result, plot_path)

    print_result(result, args.fine_theta_grid, plot_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
