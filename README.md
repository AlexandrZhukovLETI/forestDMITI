# Forest Rectangle Search

Поиск короткой ломаной для задачи выхода из прямоугольного леса.

---

## Задача программы

Поиск короткой ломаной для задачи «Как выйти из леса?» в прямоугольнике F_r=[0,r]x[0,1].

## Идея программы:
  1) Ломаная задается n равными звеньями. Сначала считаем длину звена равной 1.
  2) Для формы ломаной вычисляем, во сколько раз ее надо масштабировать,
     чтобы она гарантированно не помещалась в прямоугольник при любом повороте.
  3) Минимизируем L = n * ell по углам звеньев.

## Установка

### Windows (PowerShell)

Клонировать проект:

```powershell
git clone https://github.com/AlexandrZhukovLETI/forestDMITI.git

cd forestDMITI
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
git clone https://github.com/AlexandrZhukovLETI/forestDMITI.git

cd forestDMITI
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
python3 main.py
```

---

## Примеры запуска

Флаг --plot являеться необязательным. Его установка позваляет сохранить график иследования в формате png.

Основной шаблон команды:

```bash
python main.py --n <число_звеньев> --r <отношение_сторон> --mode <режим> --generations <число_поколений> --theta-grid <число_углов> --refine-theta --plot
```

Обычный поиск:

```bash
python main.py --n 6 --r 1.5 --mode convex --plot
```

Более точный поиск:

```bash
python main.py --n 8 --theta-grid 1440 --refine-theta
```

Поиск среди большого количества параметров:

```bash
python main.py --scan --n-values 2,4,6,8 --r-values 1,1.25,1.5,1.75,2,2.25,2.5,3 --generations 100 --theta-grid 1440 --refine-theta --plot
```