import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# Отключает просмотр рисунков в ноутбуке
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import MaxNLocator

from opticalcoating.design_class import Design
from opticalcoating.statistics_info import *
from opticalcoating.simulation_info import *
from opticalcoating.correlation_plots import *

from pathlib import Path
from os import listdir
import csv


def appr_curve(x, y, x2, y2, fname):
    font_properties = {'size': 38,
                       'family': 'Times New Roman'}
    rc('font', **font_properties)

    ax = plt.figure(figsize=(16, 9)).gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.bar(x, y)
    plt.plot(x2, y2, 'm')

    plt.xlim(1 - 0.5, len(x) + 0.5)
    plt.ylim(0., 10.5)

    plt.xlabel('Number of layer')
    plt.ylabel('RMS error, nm')

    plt.savefig(fname)
    plt.close()


def app_func(x, C, alpha):
    return C * np.exp(alpha * (x - 1))


def processing_stat(num_stat, dir_name):
    fname_stat = 'Data_to_process/' + dir_name + '/Statistics/Statistic' + str(num_stat).zfill(3) + '.json'
    with open(fname_stat, 'r') as file:
        info = json.load(file)

    n_stats = np.sum(
        [(1 if fname[-5:] == '.json' else 0) for fname in listdir('Data_to_process/' + dir_name + '/Statistics')])

    errors = pd.DataFrame(info['error list'])
    M, N = errors.shape
    errors_rms = pd.Series([np.linalg.norm(errors.iloc[:, j]) / np.sqrt(M) for j in range(N)])

    y_value = errors_rms
    x_value = pd.Series(np.array(range(1, len(y_value) + 1), dtype=float))

    args, _ = curve_fit(app_func, x_value, y_value, bounds=((np.min(y_value), -np.inf), (np.max(y_value), np.inf)))
    C = args[0]
    alpha = args[1]

    # Построим и сохраним картинку
    Path('Data_to_process/' + dir_name + '/Pictures').mkdir(exist_ok=True)
    fname_pic = 'Data_to_process/' + dir_name + '/Pictures/Picture' + str(num_stat).zfill(3) + '.png'
    appr_curve(x_value, y_value, x_value, app_func(x_value, C, alpha), fname=fname_pic)

    des = Design(name=info['design'])
    d_th = des.d[1:]
    delta = 100 * np.linalg.norm(errors_rms / d_th) / np.sqrt(N)

    fname_c = 'Data_to_process/' + dir_name + '/c_values/c_value' + str(num_stat).zfill(3) + '.json'
    with open(fname_c, 'r') as file:
        c_info = json.load(file)

    return alpha, C, delta, c_info['c_value']


def processing_all_stat(dir_name):
    n_stats = np.sum([(1 if fname[-5:] == '.json' else 0)
                      for fname in listdir('Data_to_process/' + dir_name + '/Statistics')])

    alpha_value = np.empty(n_stats, dtype=float)
    C_value = np.empty(n_stats, dtype=float)
    delta_value = np.empty(n_stats, dtype=float)
    c_value = np.empty(n_stats, dtype=float)

    fname_b = 'Data_to_process/' + dir_name + '/betas.txt'
    with open(fname_b, newline='') as csvfile:
        betas_txt = csv.reader(csvfile)
        beta_value = np.array([float(''.join(row)) for row in betas_txt])

    for num_stat in range(1, n_stats + 1):
        (alpha_value[num_stat - 1],
         C_value[num_stat - 1],
         delta_value[num_stat - 1],
         c_value[num_stat - 1]) = processing_stat(num_stat, dir_name)

    return alpha_value, C_value, delta_value, c_value, beta_value


def process_statistics(dir_name):
    alpha_value, C_value, delta_value, c_value, beta_value = processing_all_stat(dir_name)
    df = pd.DataFrame({'C': C_value,
                       'alpha': alpha_value,
                       'delta': delta_value,
                       'c': c_value,
                       'beta': beta_value})
    df.transpose().to_excel('Data_to_process/' + dir_name + '/' + dir_name + '_processed.xlsx')
