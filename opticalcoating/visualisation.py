import plotly.graph_objects as go
import pandas as pd
import seaborn as sns
import numpy as np
import json
from matplotlib import rc
import matplotlib.pyplot as plt
from .units import Wave
from .calc_flux import calc_flux


rc('font', size=30, family='Times New Roman')

# --------------------------
def bar(heights, labels=None, ymax=None, color=None, special_layers=(), fname=None, **kwargs):
    if 'figsize' not in kwargs: kwargs['figsize'] = (16, 9)
    fig = plt.figure(figsize=kwargs['figsize'])
    ax = fig.add_subplot()

    N = len(heights)
    if color is None:
        color = []
        for i in range(N):
            if i + 1 in special_layers: color.append('m')
            else: color.append('r' if i % 2 else 'b')

    ax.bar(range(1, N + 1), heights, color=color)
    plt.xlim(1 - 0.5, N + 0.5)
    plt.ylim(0., 1.05 * max(heights) if ymax is None else ymax)
    if labels is not None:
        plt.xlabel(labels['x'])
        plt.ylabel(labels['y'])
    if fname is not None: plt.savefig(fname)


def thickness_bar(des, lang='en', pic_ext=None, **kwargs):
    labels = {'ru': {'x': 'Номер слоя',
                     'y': 'Физическая толщина, нм'},
              'en': {'x': 'Layer number',
                     'y': 'Thickness, nm'}}
    kwargs['heights'] = des.d[1:des.N + 1]
    kwargs['labels'] = labels[lang]
    if pic_ext: kwargs['fname'] = f'thicknesses_{des.name}.{pic_ext}'
    bar(**kwargs)


def rms_bar(stat, ymax=None, lang='en', pic_ext=None, **kwargs):
    labels = {'ru': {'x': 'Номер слоя',
                     'y': 'Среднеквадратичная ошибка на слое, нм'},
              'en': {'x': 'Layer number',
                     'y': 'RMS error, nm'}}
    kwargs['heights'] = stat.error_rms()
    kwargs['labels'] = labels[lang]
    if ymax: kwargs['ymax'] = ymax
    if pic_ext:
        wv = int(min(stat.waves[1:]))
        algs = ''.join(tuple(set(stat.term_algs[1:])))
        meas_sigma = str(100 * max(stat.meas_sigmas[1:])).replace('.', '')
        kwargs['fname'] = f'rms_{stat.des_name}_{algs}_T{meas_sigma}_{wv}nm.{pic_ext}'
    bar(**kwargs)

# --------------------------
def plot(wv_range, values_list, labels=None, fname=None, **kwargs):
    if 'figsize' not in kwargs: kwargs['figsize'] = (16, 9)
    fig = plt.figure(figsize=kwargs['figsize'])
    ax = fig.add_subplot()

    for y_values in values_list:
        ax.plot(wv_range, y_values)

    plt.xlim(wv_range[0], wv_range[-1])
    if 'ylim' in kwargs: plt.ylim(kwargs['ylim'])
    elif '%' in labels['y']: plt.ylim(0., 100.)

    if labels is not None:
        plt.xlabel(labels['x'])
        plt.ylabel(labels['y'])
    if fname is not None: plt.savefig(fname)


def spectral_plot(designs, *, q_TR='T', wv_bnd=None, q_subs=True, lang='en', pic_ext=None, step=0.5,  **kwargs):
    labels = {'ru': {'x': 'Длина волны, нм',
                     'y': f'{q_TR}, %'},
              'en': {'x': 'Wavelength, nm',
                     'y': f'{q_TR}, %'}}
    kwargs['labels'] = labels[lang]

    if wv_bnd is None: wv_bnd = (380, 760)
    # for des in designs:
    #     left  = max(des.wv_bnd(j)[0] for j in range(des.N + 1))
    #     left = 380 if left < 0.1 else left
    #     right = min(des.wv_bnd(j)[1] for j in range(des.N + 1))
    #     right = 760 if right > 10_000 else right
    #     wv_bnd = (max(left, wv_bnd[0]), min(right, wv_bnd[1]))
    kwargs['wv_range'] = np.linspace(*wv_bnd, num=int((1/step)*(wv_bnd[1]-wv_bnd[0])) + 1)
    if 'both' in kwargs.get('polarisation', 'S'):
        waves_S = list(map(lambda wv: Wave(wv, 'S', kwargs.get('angle', 0)),
                         kwargs['wv_range']))
        waves_P = list(map(lambda wv: Wave(wv, 'P', kwargs.get('angle', 0)),
                         kwargs['wv_range']))
        for des in designs:
            kwargs.setdefault('values_list', []).append(
                np.vectorize(lambda x: calc_flux(des, x, q_subs=q_subs, q_percent=True, q_TR=q_TR))(waves_S))
        for des in designs:
            kwargs.setdefault('values_list', []).append(
                np.vectorize(lambda x: calc_flux(des, x, q_subs=q_subs, q_percent=True, q_TR=q_TR))(waves_P))
    else:
        waves = list(map(lambda wv: Wave(wv, kwargs.get('polarisation', 'S'), kwargs.get('angle', 0)),
                         kwargs['wv_range']))
        for des in designs:
            kwargs.setdefault('values_list', []).append(
                np.vectorize(lambda x: calc_flux(des, x, q_subs=q_subs, q_percent=True, q_TR=q_TR))(waves))

    if 'statistic_num' in kwargs and pic_ext:
        kwargs['fname'] = f'comp_{"best" if kwargs["is_best"] else "worst"}_{kwargs["statistic_num"]}.{pic_ext}'
    elif pic_ext: kwargs['fname'] = f'spectral_plot_{des.name}.{pic_ext}'
    plot(**kwargs)


def dispersion_plot(des, layer_num, wv_bnd=None, lang='en', pic_ext=None, **kwargs):
    labels = {'ru': {'x': 'Длина волны, нм',
                     'y': 'n'},
              'en': {'x': 'Wavelength, nm',
                     'y': 'n'}}
    kwargs['labels'] = labels[lang]

    if wv_bnd is None: wv_bnd = (0.1, 10_000)
    left  = max(des.wv_bnd(j)[0] for j in range(des.N + 1))
    left = 380 if left < 0.1 else left
    right = min(des.wv_bnd(j)[1] for j in range(des.N + 1))
    right = 760 if right > 10_000 else right
    wv_bnd = (max(left, wv_bnd[0]), min(right, wv_bnd[1]))

    kwargs['wv_range'] = np.linspace(*wv_bnd, num=int(2*(wv_bnd[1]-wv_bnd[0])))
    waves = list(map(Wave, kwargs['wv_range']))
    kwargs['values_list'] = (np.vectorize(lambda x: des.n(layer_num, x))(waves),)
    if pic_ext: kwargs['fname'] = f'dispersion_{des.name}_{des.info["layers"][layer_num]}.{pic_ext}'

    plot(**kwargs)

# --------------------------
def c_hist(stat, *, xmax=None, lang='en', pic_ext=None):
    data = pd.Series(stat.c_array)

    if xmax is None: xmax = data.max()
    else: data = data[data < xmax]
    xmin = min(data)

    plt.figure(figsize=(16, 9))
    sns.histplot(data=data, bins=40)

    plt.xlim(0 if xmin > 0 else 1.05 * xmin, 1.05 * xmax)

    labels = {'ru': {'x': 'Значение c',
                     'y': 'Число симуляций'},
              'en': {'x': 'c',
                     'y': 'Amount of simulations'}}

    plt.xlabel(labels[lang]['x'])
    plt.ylabel(labels[lang]['y'])

    if pic_ext:
        fname = f'c_hist_{stat.statistic_num}.{pic_ext}'
        plt.savefig(fname)


def thickness_error_box_plot(*, num, title='', y_range=None, special_layers=[]):
    """Box Plot visualisation of thickness errors
    obtained during simulations
    :param num: a number of the file in the folder "Statistics"
    :param title: a title of a plot
    :param y_range: a list of boundaries of thickness error axis in nm: [min, max]
    :param special_layers: a list of highlighted layers"""

    fname = 'Statistics/Statistic' + str(num).zfill(3) + '.json'
    with open(fname, 'r') as file:
        stat = json.load(file)
    N = len(stat['error list'][0])
    df_stat= pd.DataFrame(data=stat['error list'], columns=list(range(1, N + 1)))

    fig = go.Figure()

    for j in range(1, N + 1):
        if j in special_layers:
            mcol = 'purple'
        elif j % 2:
            mcol = 'blue'
        else:
            mcol = 'red'

        fig.add_trace(go.Box(y=df_stat[j],
                             name=j,
                             marker_color=mcol))

    fig.update_layout(
        title=title,
        xaxis_title="Layer number",
        yaxis_title="Thickness error, nm")

    if y_range is not None:
        fig.update_yaxes(range=y_range)
    fig.update_yaxes(minor=dict(showgrid=True))
    return fig

def thickness_error_violin_plot(*, num, title='', y_range=None, special_layers=[]):
    """Violin Plot visualisation of thickness errors
    obtained during simulations
    :param num: a number of the file in the folder "Statistic"
    :param title: a title of a plot
    :param y_range: a list of boundaries of thickness error axis in nm: [min, max]
    :param special_layers: a list of highlighted layers"""

    fname = 'Statistics/Statistic' + str(num).zfill(3) + '.json'
    with open(fname, 'r') as file:
        stat = json.load(file)
    N = len(stat['error list'][0])
    df_stat= pd.DataFrame(data=stat['error list'], columns=list(range(1, N + 1)))

    fig = go.Figure()

    for j in range(1, N + 1):
        if j in special_layers:
            mcol = 'purple'
        elif j % 2:
            mcol = 'blue'
        else:
            mcol = 'red'

        fig.add_trace(go.Violin(y=df_stat[j],
                                name=j,
                                marker_color=mcol))

    fig.update_layout(
        title=title,
        xaxis_title="Layer number",
        yaxis_title="Thickness error, nm")

    if y_range is not None:
        fig.update_yaxes(range=y_range)
    fig.update_yaxes(minor=dict(showgrid=True))
    return fig

