import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from .calc_flux import Wave, calc_flux


def monitoring_curve_plot(des, waves, control_wv=1000, lang='en', **kwargs):
    # kwargs - аргументы для calc_flux
    grid = {j: np.linspace(0, des.d_th[j], (int(des.d_th[j] // 2.0) + 1)) for j in range(1, des.N + 1)}

    def thicknesses(j, i):
        return [des.d_th[j_] for j_ in range(j)] + [grid[j][i]] + (des.N - j) * [0.]

    def total_opt_thick(j, i):
        return sum(
            n * d for n, d in zip([des.n(j, Wave(control_wv)) for j in range(1, des.N + 1)], thicknesses(j, i)[1:]))

    X = {j: [total_opt_thick(j, i) for i in range(len(grid[j]))] for j in range(1, des.N + 1)}
    Y = {}
    for j in range(1, des.N + 1):
        for i in range(len(grid[j])):
            des.d = thicknesses(j, i)
            kwargs.setdefault('q_percent', True)
            Y.setdefault(j, []).append(calc_flux(des, waves[j], **kwargs))

    rc('font', size=30, family='Times New Roman')
    plt.figure(figsize=(16, 9))

    for j in range(1, des.N + 1):
        color_plot = 'red' if not j % 2 else 'blue'
        plt.plot(X[j], Y[j], color=color_plot)

    plt.xlim(min(X[min(X)]), max(X[max(X)]))
    plt.ylim(0., 100.)
    if lang in 'ru':
        plt.xlabel('Оптическая толщина покрытия, нм')
    elif lang in 'en':
        plt.xlabel('Optical thickness, nm')
    plt.ylabel(kwargs.get('q_TR', 'R') + ', %')