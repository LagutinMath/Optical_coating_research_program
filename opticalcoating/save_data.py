from os import listdir
from importlib.resources import files

def find_file_name(obj_name: str, ext='.json'):
    """Ищет и возвращает свободное имя в папке *имя объекта*s
    Пример: obj_name = 'Simulation', первое возвращенное имя будет:
    'Simulations/Simulation001.json'"""
    names = listdir(files(f'opticalcoating.resources.{obj_name}s'))
    nums = {int(name[len(obj_name):len(obj_name) + 3]) for name in names
            if name[:len(obj_name)] in obj_name}
    num = str(min(set(range(1, max(nums, default=0) + 2)) - nums)).zfill(3)
    return files(f'opticalcoating.resources.{obj_name}s').joinpath(f'{obj_name}{num}{ext}')