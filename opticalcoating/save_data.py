from os import listdir

def find_file_name(obj_name: str, ext='.json'):
    """Ищет и возвращает свободное имя в папке *имя объекта*s
    Пример: obj_name = 'Simulation', первое возвращенное имя будет:
    'Simulations/Simulation001.json'"""
    dir_file_names = listdir(obj_name + 's/')
    indx = 1
    file_name = obj_name + 's/' + obj_name + str(indx).zfill(3) + ext

    while indx < 1000:
        file_name = obj_name + 's/' + obj_name + str(indx).zfill(3) + ext
        is_name_in_dir = False
        for name in dir_file_names:
            is_name_in_dir = ((obj_name + 's/' + name) == file_name)
            if is_name_in_dir:
                break
        if is_name_in_dir:
            indx += 1
            continue
        else:
            break

    return file_name