import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
from .statistics_info import StatInfo
import numpy as np


def std_values(errors):
    """Roots of eig values (standart deviations) of error matrix"""
    # M - number of simulations, m - number of layers
    M, m = errors.shape
    mu = (1 / M) * errors.T @ errors
    vals, _ = np.linalg.eig(mu)
    return sorted(np.sqrt(vals), reverse=True)


def sigmas_plot(statistic_num, *, ymax=None):
    font_properties = {'size': 22,
                       'family': 'Times New Roman'}
    rc('font', **font_properties)

    error_list = StatInfo.load(statistic_num)
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
