import plotly.graph_objects as go
import pandas as pd
import json


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

