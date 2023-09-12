from os import listdir
from importlib.resources import files
import pandas as pd

from .statistics_info import StatInfo
from .c_values import ProcessedStatistics


def stat_nums():
    names = listdir(files(f'opticalcoating.resources.Statistics'))
    return {int(name[9:12]) for name in names if name.startswith('Statistic')}


def process_statistics():
    C_values, alpha_values, delta_values, c_values, beta_values = [], [], [], [], []
    for stat_num in stat_nums():
        stat = StatInfo.load(stat_num)
        pr_stat = ProcessedStatistics.load(stat_num)
        C_values.append(pr_stat.exp_appr_coef['C'])
        alpha_values.append(pr_stat.exp_appr_coef['alpha'])
        delta_values.append(stat.delta())
        c_values.append(pr_stat.c_value)
        beta_values.append(pr_stat.beta)

    df = pd.DataFrame({'C': C_values,
                       'alpha': alpha_values,
                       'delta': delta_values,
                       'c': c_values,
                       'beta': beta_values},
                      index=stat_nums())
    fname = files('opticalcoating.resources.c_values').joinpath('processed.xlsx')
    df.transpose().to_excel(fname)
