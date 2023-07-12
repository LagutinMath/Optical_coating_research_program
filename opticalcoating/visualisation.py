import plotly.graph_objects as go
import pandas as pd
import seaborn as sns
import json
from matplotlib import rc
import matplotlib.pyplot as plt


rc('font', size=22, family='Times New Roman')

# --------------------------
def bar(heights, labels=None, ymax=None, color=None, special_layers=(), fname=None):
    fig = plt.figure('Thicknesses', figsize=(16, 9))
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
    kwargs['heights'] = des.d[1:des.N]
    kwargs['labels'] = labels[lang]
    if pic_ext: kwargs['fname'] = f'thicknesses_{des.name}.{pic_ext}'
    bar(**kwargs)


def rms_bar(stat, ymax=None, lang='en', pic_ext=None, **kwargs):
    labels = {'ru': {'x': 'Номер слоя',
                     'y': 'Среднеквадратичная ошибка на слое, нм'},
              'en': {'x': 'Layer number',
                     'y': 'Root mean square, nm'}}
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

def c_hist(stat, *, xmax=None, lang='en'):
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

