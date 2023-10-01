import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from .calc_flux import Wave, calc_flux

class MonochromStrategyData():
    def __init__(self, des, waves, X, Y, plot_params, control_wv=1000):
        self.des = des
        self.waves = waves
        self.control_wv = control_wv
        self._grid = None
        self.X = X
        self.Y = Y
        self.plot_params = plot_params

    @classmethod
    def brute_force_calc(cls, des, waves, step=2.0, control_wv=1000, **plot_params):
        grid = {j: np.linspace(0, des.d_th[j], (int(des.d_th[j] // step) + 1)) for j in range(1, des.N + 1)}

        def thicknesses(j, i):
            return [des.d_th[j_] for j_ in range(j)] + [grid[j][i]] + (des.N - j) * [0.]

        def total_opt_thick(j, i):
            return sum(n * d for n, d in zip([des.n(j, Wave(control_wv)) for j in range(1, des.N + 1)],
                                             thicknesses(j, i)[1:]))

        X = {j: [total_opt_thick(j, i) for i in range(len(grid[j]))] for j in range(1, des.N + 1)}
        Y = {}
        for j in range(1, des.N + 1):
            for i in range(len(grid[j])):
                des.d = thicknesses(j, i)
                plot_params.setdefault('q_percent', True)
                Y.setdefault(j, []).append(calc_flux(des, waves[j], **plot_params))

        return cls(des, waves, X, Y, plot_params, control_wv)


    def monitoring_curve_plot(self, lang='en'):
        rc('font', size=30, family='Times New Roman')
        plt.figure(figsize=(16, 9))

        for j in range(1, self.des.N + 1):
            color_plot = 'red' if not j % 2 else 'blue'
            plt.plot(self.X[j], self.Y[j], color=color_plot)

        plt.xlim(min(self.X[min(self.X)]),
                 max(self.X[max(self.X)]))
        plt.ylim(0., 100.)
        if lang in 'ru':
            plt.xlabel('Оптическая толщина покрытия, нм')
        elif lang in 'en':
            plt.xlabel('Optical thickness, nm')
        plt.ylabel(self.plot_params.get('q_TR', 'R') + ', %')

    @staticmethod
    def monitoring_curve_plots(data, lang='en', layers=None):
        rc('font', size=30, family='Times New Roman')
        fig = plt.figure(figsize=(16, 9))
        ax = fig.add_subplot()

        if layers is None: layers = list(range(1, data[0].des.N + 1))

        for msd in data:
            for j in layers:
                color_plot = 'red' if not j % 2 else 'blue'
                ax.plot(msd.X[j], msd.Y[j], color=color_plot)

        plt.xlim(min(data[0].X[min(layers)]),
                 max(data[0].X[max(layers)]))
        plt.ylim(0., 100.)
        if lang in 'ru':
            plt.xlabel('Оптическая толщина покрытия, нм')
        elif lang in 'en':
            plt.xlabel('Optical thickness, nm')
        plt.ylabel(msd.plot_params.get('q_TR', 'R') + ', %')