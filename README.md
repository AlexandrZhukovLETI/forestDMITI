# Forest Rectangle Search

Поиск короткой ломаной для задачи выхода из прямоугольного леса.

---

## Установка

### Windows (PowerShell)

Клонировать проект:

```powershell
git clone https://github.com/ТВОЙ_НИК/forest-rectangle-search.git

cd forest-rectangle-search
```

Создать виртуальное окружение:

```powershell
python -m venv .venv
```

Активировать:

```powershell
.venv\Scripts\activate
```

Установить зависимости:

```powershell
pip install -r requirements.txt
```

Запуск:

```powershell
python DMTI.py
```

---

### Linux

Клонировать проект:

```bash
git clone https://github.com/ТВОЙ_НИК/forest-rectangle-search.git

cd forest-rectangle-search
```

Создать виртуальное окружение:

```bash
python3 -m venv .venv
```

Активировать:

```bash
source .venv/bin/activate
```

Установить зависимости:

```bash
pip install -r requirements.txt
```

Запуск:

```bash
python3 DMTI.py
```

---

## Примеры запуска

Основной шаблон команды:

```bash
python DMTI.py --n <число_звеньев> --r <отношение_сторон> --mode <режим> --generations <число_поколений> --theta-grid <число_углов> --refine-theta --plot
```

Обычный поиск:

```bash
python DMTI.py --n 6 --r 1.5 --mode convex --plot
```

Более точный поиск:

```bash
python DMTI.py --n 8 --theta-grid 1440 --refine-theta
```

Поиск среди большого количества параметров:

```bash
python forest_rectangle_search.py --scan --n-values 2,4,6,8 --r-values 1,1.25,1.5,1.75,2,2.25,2.5,3 --generations 100 --theta-grid 1440 --refine-theta --plot
```