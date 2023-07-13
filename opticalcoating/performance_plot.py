import opticalcoating.visualisation as vis
from .statistics_info import StatInfo
from .c_values import ProcessedStatistics
from .design_class import Design
import numpy as np


def performance_plot(statistic_num, amount=3, is_best=True, **kwargs):
    kwargs['statistic_num'] = statistic_num
    kwargs['is_best'] = is_best
    stat = StatInfo.load(statistic_num)
    errors = stat.error_list
    sorted_sim_indx = np.argsort(ProcessedStatistics.load(statistic_num).c_array)
    sorted_sim_indx = sorted_sim_indx[:amount] if is_best else sorted_sim_indx[-amount:]
    d_th = Design(name=stat.des_name).d

    kwargs['designs'] =  [Design(name=stat.des_name) for _ in range(amount)]
    for des, sim_num in zip(kwargs['designs'], sorted_sim_indx):
        des.d = [d_th[0]] + [d_th[i] + errors[sim_num][i - 1] for i in range(1, des.N + 1)]

    vis.spectral_plot(**kwargs)