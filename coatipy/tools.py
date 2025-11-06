from os import listdir
from importlib.resources import files
import functools, time

def timer(iters=1):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            total = 0
            for i in range(iters):
                start = time.perf_counter()
                value = func(*args, **kwargs)
                end = time.perf_counter()
                total += end - start
            print(f'Среднее время выполнения {func.__name__}: {total/iters:.2} сек.')
            return value
        return wrapper
    return decorator

def find_file_name(obj_name: str, ext='.json'):
    """Ищет и возвращает свободное имя в папке *имя объекта*s
    Пример: obj_name = 'Simulation', первое возвращенное имя будет:
    'Simulations/Simulation001.json'"""
    names = listdir(files(f'opticalcoating.resources.{obj_name}s'))
    nums = {int(name[len(obj_name):len(obj_name) + 3]) for name in names
            if name[:len(obj_name)] in obj_name}
    num = str(min(set(range(1, max(nums, default=0) + 2)) - nums)).zfill(3)
    return files(f'opticalcoating.resources.{obj_name}s').joinpath(f'{obj_name}{num}{ext}')