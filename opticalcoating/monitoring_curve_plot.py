import matplotlib.pyplot as plt
from matplotlib import rc
from .calc_flux import Wave
from .deposition_simulation import simulation, SetUpParameters


def monitoring_curve_plot(*, des, waves, q_TR='R', backside=False, q_subs=True, control_wv=1000, lang='en'):
    term_algs = (des.N + 1) * ['Elimination']
    # стратегия waves q_TR, backside, (?) witness_layers
    set_up_pars = SetUpParameters(N=des.N, waves=waves, q_TR=q_TR, backside=backside)
    sim_res = simulation(des, term_algs, set_up_pars, rnd_seed=None, q_subs=q_subs)

    font_properties = {'size': 22,
                       'family': 'Times New Roman'}
    rc('font', **font_properties)

    # несквозной список моментов времени (каждый слой начинает напыляться в момент времени t = 0)
    layer_time = [[t - time_list[0] for t in time_list] for time_list in sim_res.time_list]

    X = []
    Y = []
    plt.figure(figsize=(16, 9))

    # полная оптическая толщина
    h = 0.
    for j in range(1, des.N + 1):
        # скорость r = set_up_pars.rates[j] (0.5 нм)
        r = set_up_pars.rates[j]
        n = des.n(j, Wave(control_wv))

        X.append([h + r * n * t for t in layer_time[j]])
        h += r * n * layer_time[j][-1]

        Y.append([100. * TR for TR in sim_res.flux_meas[j]])
        if j % 2 == 0:
            color_plot = 'red'
        else:
            color_plot = 'blue'
        plt.plot(X[j - 1], Y[j - 1], color=color_plot)

    plt.xlim(min(X[0]), max(X[-1]))
    plt.ylim(0., 100.)
    if lang=='ru':
        plt.xlabel(f'Оптическая толщина покрытия (контрольная длина волны = {control_wv} нм), нм')
    else:
        plt.xlabel('Optical thickness, nm')
    plt.ylabel(q_TR + ', %')