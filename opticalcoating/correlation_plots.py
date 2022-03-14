from opticalcoating.correlation import std_values
import pandas as pd
import matplotlib.pyplot as plt
from opticalcoating.statistics_info import load_dict
from matplotlib import rc
from opticalcoating.save_data import find_file_name


def sigmas_plot(statistic_num, *, show=False, ymax=None):
    font_properties = {'size': 22,
                       'family': 'Times New Roman'}
    rc('font', **font_properties)

    error_list = load_dict(statistic_num)
    lolik = pd.DataFrame(error_list['error list'])

    stds = std_values(lolik)

    plt.figure(figsize=(16,9))
    plt.bar(x=range(1, len(stds) + 1), height=stds)
    plt.xlim(1 - 0.5, len(stds) + 0.5)
    if ymax is None:
        ymax = max(stds)
    plt.ylim(0., 1.05 * ymax)
    plt.xlabel('Номер собственного вектора')
    plt.ylabel('Стандартное отклонение ошибки, нм')
    plt.savefig(find_file_name('Picture', '.png'))

    if show:
        plt.show()

