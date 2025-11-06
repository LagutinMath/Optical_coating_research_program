from os import listdir
from importlib.resources import files
import pandas as pd

from .statistics_info import StatInfo
from .c_values import ProcessedStatistics

def essential_info(info):
    info_set = set(info)
    info_set.discard(None)
    seq = [f'{each:.3}' if isinstance(each, float) else str(each) for each in info_set]
    return ', '.join(seq)

def stat_nums():
    names = listdir(files(f'opticalcoating.resources.Statistics'))
    return {int(name[9:12]) for name in names if name.startswith('Statistic')}


def process_statistics():
    info_dct = {}
    for stat_num in stat_nums():
        stat = StatInfo.load(stat_num)
        pr_stat = ProcessedStatistics.load(stat_num)

        info_dct.setdefault('des_name', []).append(stat.des_name)
        info_dct.setdefault('N_sim', []).append(len(stat.error_list))
        info_dct.setdefault('term_algs', []).append(essential_info(stat.term_algs))
        info_dct.setdefault('width', []).append(stat.width if stat.width is not None else 0)
        info_dct.setdefault('waves', []).append(essential_info(stat.waves))
        info_dct.setdefault('rates', []).append(essential_info(stat.rates))
        info_dct.setdefault('rates_sigmas', []).append(essential_info(stat.rates_sigmas))
        info_dct.setdefault('meas_sigmas', []).append(essential_info(stat.meas_sigmas))

        info_dct.setdefault('mean_error_norm', []).append(stat.mean_error_norm)
        info_dct.setdefault('median_error_norm', []).append(stat.median_error_norm)
        info_dct.setdefault('C', []).append(pr_stat.exp_appr_coef['C'])
        info_dct.setdefault('alpha', []).append(pr_stat.exp_appr_coef['alpha'])
        info_dct.setdefault('delta', []).append(stat.delta())
        info_dct.setdefault('c', []).append(pr_stat.c_value)
        info_dct.setdefault('beta', []).append(pr_stat.beta)


    df = pd.DataFrame(info_dct, index=stat_nums())
    fname = files('opticalcoating.resources.c_values').joinpath('processed.xlsx')
    df.transpose().to_excel(fname)
