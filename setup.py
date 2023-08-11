from pathlib import Path
from importlib.resources import files

folders = ('c_values', 'Designs', 'Simulations', 'Spectral_values', 'Statistics', 'Targets')
for folder in folders:
    f = files('opticalcoating.resources').joinpath(folder)
    Path(f).mkdir(exist_ok=True)
    with open(f.joinpath('__init__.py'), 'w') as f: pass