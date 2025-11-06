import copy
from .deposition_simulation import MonochromStrategyInfo, DataNonloc

class MonochromaticDeposition():
    def __init__(self, des_th, set_up_pars, t_scan):
        self.des_th = copy.deepcopy(des_th)
        self.set_up_pars = set_up_pars
        self.str_info = MonochromStrategyInfo(des_th, set_up_pars)
        self.nonloc_alg = [None] + [DataNonloc(j, des_th, set_up_pars) for j in range(1, des_th.N + 1)]
        self.t_scan = t_scan
        self.dt = 0.


    def update(self, j, flux_meas, t_scan=None):
        """Обновление данных объекта
        :param flux_meas: измеренный на слое j коэф. проп./отр в абс. величинах T/R in [0., 1.]
        :param dt: сколько времени прошло с начала напыления слоя"""
        if t_scan is None:
            self.dt += self.t_scan
        else:
            self.dt += t_scan
        self.j = j
        self.nonloc_alg[j].refresh(self.dt, flux_meas, self.set_up_pars.q_TR[j])


    def term_time_predict(self, term_algs='Elimination'):
        if term_algs == 'Elimination':
            term_flux_lvl = self.nonloc_alg[self.j].flux(self.des_th.d[self.j], self.set_up_pars.q_TR[self.j])
        elif term_algs == 'Quasiswing':
            q_turn = self.str_info.q_prev_extr[self.j]
            if q_turn == 'max':
                turn = self.nonloc_alg[self.j].coef.flux_max(q_TR=self.set_up_pars.q_TR[self.j])
            else:
                turn = self.nonloc_alg[self.j].coef.flux_min(q_TR=self.set_up_pars.q_TR[self.j])
            ampl = self.nonloc_alg[self.j].coef.ampl()
            term_flux_lvl = turn - (self.str_info.prev_extr[self.j]
                                    - self.str_info.term[self.j]) * (ampl / self.str_info.ampl[self.j])
        elif term_algs[self.j] == 'None':
            term_flux_lvl = self.str_info.term[self.j]
        else:
            raise NameError(f'Uncorrect termination algoritm for layer {self.j}')

        t_term = self.nonloc_alg[self.j].calc_t(d_th=self.des_th.d[self.j],
                                                lvl=term_flux_lvl, q_TR=self.set_up_pars.q_TR[self.j])
        return t_term
