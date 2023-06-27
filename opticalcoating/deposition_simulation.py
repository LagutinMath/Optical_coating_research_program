import os.path
import numpy as np
import math
from opticalcoating.design_class import *
from opticalcoating.simulation_info import SimInfo
import copy
from matplotlib.animation import FuncAnimation
from datetime import datetime
from numpy import pi


class SetUpParameters:
    """Класс содержащий в себе описание
    1. Параметров установки заданных пользователем, таких как:
        - длина волны мониторинга,
        - наблюдается коэф.отражения или пропускания
        - промежуток времени через который проводятся измерения
    2. Информацию о несовершенстве установки"""

    def __init__(self, *, N, waves=510, q_TR='R', tau=2., rates=0.5, mono_width=0., meas_sigmas=None, meas_syst=0.,
                 rates_sigmas=None, rates_syst=0., delay_time=0., delay_time_sigma=0., r_index_syst=None,
                 r_index_sigmas=None, backside=False):
        self.N = N

        if isinstance(waves, (int, float)):
            self.waves = [Wave(waves) for all in range(self.N + 1)]
        else:
            self.waves = waves

        if q_TR == 'T' or q_TR == 'R':
            self.q_TR = [q_TR for iter in range(self.N + 1)]
        else:
            self.q_TR = q_TR

        self.tau = tau
        self.backside = backside

        self.rates = autofill(rates, N=self.N)
        self.mono_width = mono_width
        self.meas_sigmas = autofill(meas_sigmas, N=self.N)
        self.meas_syst = meas_syst
        self.rates_sigmas = autofill(rates_sigmas, N=self.N)
        self.rates_syst = rates_syst
        self.delay_time = delay_time
        self.delay_time_sigma = delay_time_sigma
        self.r_index_syst = autofill(r_index_syst, N=self.N)
        self.r_index_sigmas = autofill(r_index_sigmas, N=self.N)
        # N --- кол-во слоёв
        # в массивах индекс соответствует номеру слоя
        # индекс 0 незадействован для rates и пр.


    def strategy_spreadsheet_plot(self):
        """Построение кривой монохроматического мониторинга
        """
        pass

def autofill(x, *, N):
    if x is None:
        return (N + 1) * [0.]
    else:
        if isinstance(x, (int, float)):
            return (N + 1) * [x]
        else:
            return x

def num_step_estimation(des: Design, set_up_pars: SetUpParameters):
    """Вычисление оценки кол-ва шагов в предстоящей симуляции"""
    r = set_up_pars.rates
    tau = set_up_pars.tau
    total_time = 0.0
    for j in range(1, des.N + 1):
        total_time += des.d[j] / r[j]
    return 10 * math.ceil(total_time / tau)

def num_step_layer_estimation(des: Design, set_up_pars: SetUpParameters):
    """Возвращает лист нужный для генерации пустых списков для предстоящей симуляции
    первое число кол-во слоёв (кол-во строк)
    второе число кол-во шагов на одном слое (кол-во столбцов)
    кол-во шагов возвращает с 10-ти кратным запасом"""
    r = set_up_pars.rates
    tau = set_up_pars.tau
    max_step = 0
    for j in range(1, des.N + 1):
        max_step = max(max_step, des.d[j] / r[j])
    max_step = math.ceil(max_step / tau)
    return des.N + 1, 10 * max_step

def norm_3sigma_rnd(rng, *, mean=0.0, sigma=0.0):
    """Функция возвращает нормально распределенную случайную величину,
    с ненулевым мат.ожиданием, которая модифицирована так, чтобы
    значение не выходило за границы интервала [mean-3*sigma, mean+3*sigma]"""
    rnd_val = rng.normal(mean, sigma)
    if rnd_val > mean + 3.0 * sigma:
        rnd_val = mean + 3.0 * sigma
    elif rnd_val < mean - 3.0 * sigma:
        rnd_val = mean - 3.0 * sigma
    return rnd_val

def my_atan(x: float, y: float):
    if x > 0:
        if y >= 0:
            return math.atan(y / x)
        else:
            return 2 * np.pi + math.atan(y / x)
    elif x < 0:
        return np.pi + math.atan(y / x)
    else:  # (x == 0)
        return 1.5 * np.pi if (y < 0) else 0.5 * np.pi

def num_cutoffs(theta: float, phi: float):
    """Количество точек pi * n принадлежащих (theta, theta + 2 * phi), где n --- целое"""
    prec = 0.0001
    b = theta / np.pi
    e = (theta + 2 * phi) / np.pi
    # для привдения к интервалу [0, pi]
    if b > e:
        b, e = e, b
        q_reverse = True
    else:
        q_reverse = False
    # эта часть выбрасывает точки на границах
    if abs(b - round(b)) < prec:
        b = round(b) + prec
    if abs(e - round(e)) < prec:
        e = round(e) - prec
    ans = math.floor(e) - math.ceil(b) + 1
    if q_reverse:
        ans = - ans - 1
    return ans


class NonlocCoef:
    def __init__(self, abg=None, Dthg=None):
        """Нелокальные коэф. можно определить
        либо как abg = (alpha, beta, gamma),
        либо как Dthg = (D, theta, gamma)"""
        if abg is None:
            self.abg = [1.0, 1.0, 1.0]
        else:
            self.abg = abg
        self.Dthg = Dthg

        self.calc_Dthg()

    def alpha(self):
        if self.abg is None:
            self.calc_abg()
        return self.abg[0]

    def beta(self):
        if self.abg is None:
            self.calc_abg()
        return self.abg[1]

    def gamma(self):
        if self.abg is None:
            return self.Dthg[2]
        else:
            return self.abg[2]

    def D(self):
        if self.Dthg is None:
            self.calc_Dthg()
        return self.Dthg[0]

    def theta(self):
        if self.Dthg is None:
            self.calc_Dthg()
        return self.Dthg[1]

    def calc_Dthg(self):
        # D = sqrt(alpha^2 + beta^2)
        D = math.sqrt(self.abg[0] ** 2 + self.abg[1] ** 2)
        # theta = my_atan(alpha / Delta, -beta / Delta);
        theta = my_atan(self.abg[0] / D, -self.abg[1] / D)
        self.Dthg = [D, theta, self.abg[2]]

    def calc_abg(self):
        # пока нигде не нужна
        pass

    def flux_max(self, q_TR='T', q_percent=False):
        if self.Dthg is None:
            self.calc_Dthg()
        res = (1.0 + 1.0 / (self.Dthg[2] - self.Dthg[0])) if (q_TR == 'R') else (-1.0 / (self.Dthg[2] + self.Dthg[0]))
        return 100.0 * res if q_percent else res

    def flux_min(self, q_TR='T', q_percent=False):
        if self.Dthg is None:
            self.calc_Dthg()
        res = (1.0 + 1.0 / (self.Dthg[2] + self.Dthg[0])) if (q_TR == 'R') else (-1.0 / (self.Dthg[2] - self.Dthg[0]))
        return 100.0 * res if q_percent else res

    def ampl(self):
        return abs(self.flux_max() - self.flux_min())


    def prev_extr(self, phi_act, q_TR):
        self.calc_Dthg()
        if math.floor((2 * phi_act + self.theta()) / np.pi) % 2 == 1:
            if q_TR == 'R':
                return self.flux_max(q_TR)
            elif q_TR == 'T':
                return self.flux_min(q_TR)
        else:
            if q_TR == 'R':
                return self.flux_min(q_TR)
            elif q_TR == 'T':
                return self.flux_max(q_TR)


    def q_prev_extr(self, phi_act, q_TR):
        self.calc_Dthg()
        if math.floor((2 * phi_act + self.theta()) / np.pi) % 2 == 1:
            if q_TR == 'R':
                return 'max'
            elif q_TR == 'T':
                return 'min'
        else:
            if q_TR == 'R':
                return 'min'
            elif q_TR == 'T':
                return 'max'

def R_full(R_1, R_2, d_s, xi, wvlen):
    """Полное отражение нормально падающего света от системы c двумя поверхностями,
    между которыми не происходит интерференции
    :param R_1: отражение от передней поверхности
    :param R_2: отражение от задней поверхности
    :param d_s: растояние между поверхностями
    :param xi: поглощение в среде между поверхностями"""
    sgm = math.exp(- 4 * pi * d_s * xi / wvlen)
    return R_1 + (sgm**2 * R_2 * (1 - R_1)**2)/(1 - sgm**2 * R_1 * R_2)

def T_full(T_1, T_2, d_s, xi, wvlen):
    """Полное пропускание нормально падающего света системы c двумя поверхностями,
    между которыми не происходит интерференции
    :param T_1: пропускание передней поверхности
    :param T_2: пропускание задней поверхности
    :param d_s: растояние между поверхностями
    :param xi: поглощение в среде между поверхностями"""
    sgm = math.exp(- 4 * pi * d_s * xi / wvlen)
    return sgm * T_1 * T_2 / (1 - sgm**2 * (1 - T_1) * (1 - T_2))

def theor_nonloccoef(A_re: float, A_im: float, n_j: float, n_s: float, only_front_side=False):
    """Вычисляет нелок.коэф-ты по адмиттансу (с учетом отражения от задней стороны подложки)
    :param A_re: Действительная часть адмиттанса
    :param A_im: Мнимая часть адмиттанса
    :param n_j: показатель преломления текущего слоя
    :param n_s: показатель преломления
    :return: Нелокальные коэффициенты"""
    if only_front_side:
        k_R = 0.
    else:
        k_R = ((n_s - 1.) / (n_s + 1)) ** 2
    rho = math.sqrt(A_re ** 2 + A_im ** 2)

    alpha = (rho ** 2 - n_j ** 2) * (1. - n_j ** 2) / (8. * A_re * n_j ** 2)
    beta = (1 - n_j ** 2) * A_im / (4. * n_j * A_re)
    gamma = 0.5 + (1. / (k_R - 1.)) - ((rho ** 2 + n_j ** 2) * (1. + n_j ** 2) / (8. * A_re * n_j ** 2))

    return NonlocCoef(abg=[alpha, beta, gamma])


class DataNonloc:
    def __init__(self, j, des, set_up_pars):
        self.wavelength = set_up_pars.waves[j].wavelength
        self.n = des.n(j, set_up_pars.waves[j])
        self.r = set_up_pars.rates[j]
        self.L = np.zeros(6, dtype=float)
        self.f = np.zeros(3, dtype=float)
        self.coef = NonlocCoef()
        self.total_refresh = 0
        self.exist_coef = False

    def refresh(self, dt, flux_meas, q_TR):
        self.total_refresh += 1

        d_phi = self.n * self.r * dt * 2 * np.pi / self.wavelength
        x = math.cos(2 * d_phi)
        y = math.sin(2 * d_phi)
        w = flux_meas if q_TR == 'T' else (1.0 - flux_meas)  # flux_meas have to be absolute value

        self.L[0] = self.L[0] + w ** 4 * x * x
        self.L[1] = self.L[1] + w ** 4 * x * y
        self.L[2] = self.L[2] + w ** 4 * x
        self.L[3] = self.L[3] + w ** 4 * y * y
        self.L[4] = self.L[4] + w ** 4 * y
        self.L[5] = self.L[5] + w ** 4

        self.f[0] = self.f[0] - w ** 3 * x
        self.f[1] = self.f[1] - w ** 3 * y
        self.f[2] = self.f[2] - w ** 3

        if self.total_refresh > 3:
            self.coef.abg[0] = (self.f[0] * (self.L[4] * self.L[4] - self.L[3] * self.L[5])
                                + self.f[1] * (self.L[1] * self.L[5] - self.L[2] * self.L[4])
                                + self.f[2] * (self.L[2] * self.L[3] - self.L[1] * self.L[4])) / (
                                       self.L[2] * self.L[2] * self.L[3]
                                       - 2 * self.L[1] * self.L[2] * self.L[4] + self.L[0] * self.L[4] * self.L[4]
                                       + self.L[1] * self.L[1] * self.L[5] - self.L[0] * self.L[3] * self.L[5])

            self.coef.abg[1] = (self.f[0] * (self.L[1] * self.L[5] - self.L[2] * self.L[4])
                                + self.f[1] * (self.L[2] * self.L[2] - self.L[0] * self.L[5])
                                + self.f[2] * (self.L[0] * self.L[4] - self.L[1] * self.L[2])) / (
                                       self.L[2] * self.L[2] * self.L[3]
                                       - 2 * self.L[1] * self.L[2] * self.L[4]
                                       + self.L[0] * self.L[4] * self.L[4]
                                       + self.L[1] * self.L[1] * self.L[5]
                                       - self.L[0] * self.L[3] * self.L[5])

            self.coef.abg[2] = (self.f[0] * (self.L[2] * self.L[3] - self.L[1] * self.L[4])
                                + self.f[1] * (self.L[0] * self.L[4] - self.L[1] * self.L[2])
                                + self.f[2] * (self.L[1] * self.L[1] - self.L[0] * self.L[3])) / (
                                       self.L[2] * self.L[2] * self.L[3]
                                       - 2 * self.L[1] * self.L[2] * self.L[4]
                                       + self.L[0] * self.L[4] * self.L[4]
                                       + self.L[1] * self.L[1] * self.L[5]
                                       - self.L[0] * self.L[3] * self.L[5])

            if not (math.isnan(self.coef.abg[0]) or math.isnan(self.coef.abg[1]) or math.isnan(self.coef.abg[2])):
                self.exist_coef = True
                # (!) изменение коэф. alpha beta gamma должно вызывать пересчёт всех коэф. в NonlocCoef
                self.coef.calc_Dthg()


    def flux(self, thickness, q_TR):
        """Возвращает значение flux, которое соответствует напыленной толщине thickness
        При подстановки d_j^{theor} это есть уровень остановки независимого метода контроля"""
        d_phi = self.n * thickness * 2 * np.pi / self.wavelength
        if q_TR == 'R':
            ans = 1.0 + 1.0 / (self.coef.D() * math.cos(2 * d_phi + self.coef.theta()) + self.coef.gamma())
        else:
            ans = -1.0 / (self.coef.D() * math.cos(2 * d_phi + self.coef.theta()) + self.coef.gamma())
        return ans


    def q_expected_interval(self, dt, d_th, *, theta_th):
        """Сообщает находится ли в момент времени dt
        процесс напыления в промежутке между теми экстремумами,
        где должен прекратиться процесс напыления
        :param dt: время прошедшее с начала напыления слоя
        :param d_th: толщина теоретического дизайна
        :param theta_th: оценка теоретического theta на слое"""
        phi_act = self.n * self.r * dt * 2 * np.pi / self.wavelength
        phi_th = self.n * d_th * 2 * np.pi / self.wavelength

        if num_cutoffs(theta_th, phi_act) == num_cutoffs(theta_th, phi_th):
            return True
        else:
            return False

    def prev_extr(self, d, q_TR):
        phi_act = self.n * d * 2 * np.pi / self.wavelength
        return self.coef.prev_extr(phi_act, q_TR)


    # def q_pass_term_flux_lvl(self, dt, term_flux_lvl, q_TR):
    #     flux_cur = self.flux(self.r * dt, q_TR)
    #     flux_prev_extr = self.prev_extr(dt, q_TR)
    #     return ((flux_prev_extr - term_flux_lvl) * (flux_cur - term_flux_lvl)) < 0.0


    def calc_t(self, *, d_th, lvl, q_TR):
        """Вычисление времени остановки по уровню прекращения напыления"""
        phi_th = self.n * d_th * 2 * pi / self.wavelength
        if q_TR == 'R':
            temp = (1.0 / (lvl - 1.0) - self.coef.gamma()) / self.coef.D()
        else:
            temp = (1.0 / (- lvl) - self.coef.gamma()) / self.coef.D()

        # Без коррекции уровня прекращения возможна ситуация при которой сигнал никогда не достигнет уровня lvl
        # в этом случае остановимся на экстремуме
        if temp > 1:
            temp = 1
        elif temp < - 1:
            temp = - 1

        # + solution + n * pi
        phi_p = 0.5 * (math.acos(temp) - self.coef.theta())
        # - solution + n * pi
        phi_n = - 0.5 * (math.acos(temp) + self.coef.theta())

        n_p = math.ceil((phi_th - phi_p) / pi)
        n_n = math.ceil((phi_th - phi_n) / pi)

        phi_p1 = phi_p + pi * (n_p - 1)
        phi_p2 = phi_p + pi * n_p
        phi_n1 = phi_n + pi * (n_n - 1)
        phi_n2 = phi_n + pi * n_n
        phi_lst = [phi_p1, phi_p2, phi_n1, phi_n2]
        dist_phi = [math.fabs(each_phi - phi_th) for each_phi in phi_lst]

        indx, elem = min(enumerate(dist_phi), key=lambda x: x[1])
        # Ближайщее решение
        phi = phi_lst[indx]

        return phi * self.wavelength / (self.n * self.r * 2 * pi)


    # def calc_delta_t(self, dt, term_flux_lvl, q_TR, tau=2.0):
    #     phi_act = self.n * self.r * dt * 2 * pi / self.wavelength
    #     if q_TR == 'R':
    #         temp = (1.0 / (term_flux_lvl - 1.0) - self.coef.gamma()) / self.coef.D()
    #     else:
    #         temp = (1.0 / (- term_flux_lvl) - self.coef.gamma()) / self.coef.D()
    #     # + solution + n * pi
    #     phi_p = 0.5 * (math.acos(temp) - self.coef.theta())
    #     # - solution + n * pi
    #     phi_n = - 0.5 * (math.acos(temp) + self.coef.theta())
    #
    #     n_p = math.ceil((phi_act - phi_p) / pi)
    #     n_n = math.ceil((phi_act - phi_n) / pi)
    #
    #     phi_p = phi_p + pi * n_p
    #     phi_n = phi_n + pi * n_n
    #     # Ближайщее решение большее, чем phi_act
    #     phi = min(phi_p, phi_n)
    #
    #     delta_t = (phi - phi_act) * self.wavelength / (self.n * self.r * 2 * pi)
    #
    #     if (delta_t >= 0.) and (delta_t <= tau):
    #         return delta_t
    #     else:
    #         return 0.0

def increase_admittance(A_re, A_im, wavelength, d, n):
    phi = (2. * np.pi / wavelength) * d * n
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)

    A_re, A_im = (
        A_re * n ** 2 / (
                (A_re * sin_phi) ** 2 + (n * cos_phi - A_im * sin_phi) ** 2
        ),

        n * (n * A_im * (cos_phi ** 2 - sin_phi ** 2)
             + cos_phi * sin_phi * (n ** 2 - A_re ** 2 - A_im ** 2)) / (
                (A_re * sin_phi) ** 2 + (n * cos_phi - A_im * sin_phi) ** 2
        )
    )
    return A_re, A_im


class MonochromStrategyInfo:
    def __init__(self, des_th: Design, set_up_pars: SetUpParameters, q_subs=True):
        """Содержит информацию о значениях энергетических коэф. T/R
        в начале напыления слоя, в конце напыления слоя и экстремумы,
        которые могут достигаться при достаточной толщине слоя
        для заданной стратегии монохроматического мониторинга
        в случае теоретического дизайна с однородными слоями"""
        self.start = (des_th.N + 1) * [0.0]
        self.max = (des_th.N + 1) * [0.0]
        self.min = (des_th.N + 1) * [0.0]
        self.term = (des_th.N + 1) * [0.0]
        self.ampl = (des_th.N + 1) * [0.0]
        self.prev_extr = (des_th.N + 1) * [0.0]
        self.q_prev_extr = (des_th.N + 1) * [0.0]
        self.num_tp = (des_th.N + 1) * [0.0]

        self.nonloc_coef = (des_th.N + 1) * [[]]

        for j_cur in range(1, des_th.N + 1):
            wavelength = set_up_pars.waves[j_cur].wavelength
            n_s = des_th.n(0, set_up_pars.waves[j_cur])
            if set_up_pars.backside or not q_subs:
                R_1 = ((n_s - 1.) / (n_s + 1)) ** 2
                T_1 = 1. - R_1
            A_re = n_s
            A_im = 0.
            for j in des_th.witness_layers(j_cur):
                d = des_th.d[j]
                n = des_th.n(j, set_up_pars.waves[j_cur])
                A_re, A_im = increase_admittance(A_re, A_im, wavelength, d, n)

            n_j_cur = des_th.n(j_cur, set_up_pars.waves[j_cur])
            if set_up_pars.backside or not q_subs:
                abg = theor_nonloccoef(A_re, A_im, n_j_cur, n_s, only_front_side=True)
                d_s = des_th.d[0]
                xi = des_th.xi(j_cur, wavelength)
                wv = wavelength
            else:
                abg = theor_nonloccoef(A_re, A_im, n_j_cur, n_s)
            abg.calc_Dthg()
            self.nonloc_coef[j_cur] = [abg.alpha(), abg.beta(), abg.gamma(), abg.D(), abg.theta()]

            if set_up_pars.q_TR[j_cur] == 'T':
                self.start[j_cur] = - 1. / (abg.alpha() + abg.gamma())
                self.max[j_cur] = - 1. / (abg.gamma() + abg.D())
                self.min[j_cur] = - 1. / (abg.gamma() - abg.D())
                if set_up_pars.backside:
                    self.start[j_cur] = T_full(T_1, self.start[j_cur], d_s, xi, wv)
                    self.max[j_cur] = T_full(T_1, self.max[j_cur], d_s, xi, wv)
                    self.min[j_cur] = T_full(T_1, self.min[j_cur], d_s, xi, wv)
            else:
                self.start[j_cur] = 1. + 1. / (abg.alpha() + abg.gamma())
                self.max[j_cur] = 1. + 1. / (abg.gamma() - abg.D())
                self.min[j_cur] = 1. + 1. / (abg.gamma() + abg.D())
                if set_up_pars.backside:
                    self.start[j_cur] = R_full(R_1, self.start[j_cur], d_s, xi, wv)
                    self.max[j_cur] = R_full(R_1, self.max[j_cur], d_s, xi, wv)
                    self.min[j_cur] = R_full(R_1, self.min[j_cur], d_s, xi, wv)

            self.ampl[j_cur] = self.max[j_cur] - self.min[j_cur]

            d_j_cur = des_th.d[j_cur]
            phi = (2. * np.pi / wavelength) * d_j_cur * n_j_cur

            if set_up_pars.q_TR[j_cur] == 'T':
                self.term[j_cur] = - 1. / (abg.alpha() * math.cos(2. * phi) + abg.beta() * math.sin(2. * phi) + abg.gamma())
                if set_up_pars.backside:
                    self.term[j_cur] = T_full(T_1, self.term[j_cur], d_s, xi, wv)
            else:
                self.term[j_cur] = 1. + 1. / (abg.alpha() * math.cos(2. * phi) + abg.beta() * math.sin(2. * phi) + abg.gamma())
                if set_up_pars.backside:
                    self.term[j_cur] = R_full(R_1, self.term[j_cur], d_s, xi, wv)
            self.num_tp[j_cur] = num_cutoffs(abg.theta(), phi)

            self.prev_extr[j_cur] = abg.prev_extr(phi, set_up_pars.q_TR[j_cur])
            if set_up_pars.backside:
                if set_up_pars.q_TR[j_cur] == 'T':
                    self.prev_extr[j_cur] = T_full(T_1, self.prev_extr[j_cur], d_s, xi, wv)
                else:
                    self.prev_extr[j_cur] = R_full(R_1, self.prev_extr[j_cur], d_s, xi, wv)
            self.q_prev_extr[j_cur] = abg.q_prev_extr(phi, set_up_pars.q_TR[j_cur])

def simulation(des_th, term_algs, set_up_pars, rnd_seed=None, q_subs=True, save_M=False):
    # save_M ускоряет счет (усл: без смены свидетелей, без смены длины волны)
    if rnd_seed is None:
        rnd_seed = datetime.now().time().microsecond
    str_info = MonochromStrategyInfo(des_th, set_up_pars, q_subs=q_subs)
    rng = np.random.default_rng(rnd_seed)
    # Для производительности определим и выделим заранее достаточное кол-во памяти массивам
    len_sim_list = num_step_layer_estimation(des_th, set_up_pars)
    max_steps = len_sim_list[1]

    # time_list = np.empty(len_sim_list, dtype=float)
    # flux_meas = np.empty(len_sim_list, dtype=float)
    # flux_act = np.empty(len_sim_list, dtype=float)

    time_list = [[0.0] for _ in range(des_th.N + 1)]
    d_j_act_t = [[0.0] for _ in range(des_th.N + 1)]
    flux_meas = [[0.0] for _ in range(des_th.N + 1)]
    flux_act = [[0.0] for _ in range(des_th.N + 1)]

    # Информация о дизайне в текущий для симуляции момент времени
    des_act = copy.deepcopy(des_th)
    # des_act.d = np.zeros(des_th.N + 1, dtype=float)
    des_act.d = (des_th.N + 1) * [0.0]
    # term_cond_case = np.zeros(des_th.N + 1, dtype=float)
    term_cond_case = (des_th.N + 1) * [0]

    for j in range(1, des_th.N + 1):
        # des_act.n_fix[0][j] = des_act.n(j, set_up_pars.waves[j])\
        #         * (1 + norm_3sigma_rnd(rng, mean=set_up_pars.r_index_syst[j], sigma=set_up_pars.r_index_sigmas[j]))
        # На новом слое получем теоретическое theta (можно просто оценить, но мы получим точно)
        nonloc_coef_th = MonochromStrategyInfo(des_act, set_up_pars, q_subs=q_subs).nonloc_coef[j]
        theta_th = nonloc_coef_th[-1]

        # print('j =', j)
        nonloc_alg = DataNonloc(j, des_th, set_up_pars)
        term_cond = False

        # в конце напыления шаг по времени может быть меньше tau
        # поэтому введём переменную delta_t означающую размер интервала времени для выполняемого шага
        delta_t = set_up_pars.tau
        # dt --- времени прошло с начала напыления слоя
        dt = 0.0
        expected_interval = False

        # шаг *измерения и анализа*
        # lsc = 0  # layer_step_counter
        time_list[j][0] = time_list[j - 1][-1]
        flux_act[j][0] = des_act.calc_flux(set_up_pars.waves[j], q_TR=set_up_pars.q_TR[j], layer=j,
                                           backside=set_up_pars.backside, q_subs=q_subs, save_M=save_M)
        flux_meas[j][0] = flux_act[j][-1] + norm_3sigma_rnd(rng, sigma=set_up_pars.meas_sigmas[j])
        nonloc_alg.refresh(dt, flux_meas[j][-1], set_up_pars.q_TR[j])

        # напыление слоя
        # while not term_cond:
        for lsc in range(1, max_steps):
            if term_algs[j] == 'Relative Thickness Error':
                term_cond_case[j] = 1
                des_act.increase_layer_thickness(j, des_th.d[j] * (1. + norm_3sigma_rnd(rng, sigma=0.01)))
                break

            # шаг *изменения состояния природы*
            dt += delta_t
            cur_rate = set_up_pars.rates[j] + norm_3sigma_rnd(rng, sigma=set_up_pars.rates_sigmas[j])
            des_act.increase_layer_thickness(j, delta_t * cur_rate)

            # шаг *измерения и анализа*
            time_list[j].append(time_list[j][-1] + delta_t)
            d_j_act_t[j].append(des_act.d[j])
            flux_act[j].append(des_act.calc_flux(set_up_pars.waves[j], q_TR=set_up_pars.q_TR[j], layer=j,
                                                 backside=set_up_pars.backside, q_subs=q_subs, save_M=save_M))
            flux_meas[j].append(flux_act[j][-1] + norm_3sigma_rnd(rng, sigma=set_up_pars.meas_sigmas[j]))
            nonloc_alg.refresh(dt, flux_meas[j][-1], set_up_pars.q_TR[j])

            if (not nonloc_alg.exist_coef) or (des_act.d[j] < 0.7 * des_th.d[j]):
                continue
            else:
                if term_algs[j] in ('Elimination', 'TM'):
                    term_flux_lvl = nonloc_alg.flux(des_th.d[j], set_up_pars.q_TR[j])
                elif term_algs[j] in ('Quasiswing', 'QS'):
                    # turn = nonloc_alg.prev_extr(des_th.d[j], set_up_pars.q_TR[j])
                    q_turn = str_info.q_prev_extr[j]
                    if q_turn == 'max':
                        turn = nonloc_alg.coef.flux_max(q_TR=set_up_pars.q_TR[j])
                    else:
                        turn = nonloc_alg.coef.flux_min(q_TR=set_up_pars.q_TR[j])
                    ampl = nonloc_alg.coef.ampl()
                    term_flux_lvl = turn - (str_info.prev_extr[j] - str_info.term[j]) * (ampl / str_info.ampl[j])
                # elif term_algs[j] == 'Pseudoswing':
                #     # Не готово
                #     lvl = 0.0
                elif term_algs[j] in ('None', 'LM'):
                    term_flux_lvl = str_info.term[j]
                else:
                    raise NameError(f'Uncorrect termination algoritm for layer {j}')

                t_term = nonloc_alg.calc_t(d_th=des_th.d[j], lvl=term_flux_lvl, q_TR=set_up_pars.q_TR[j])

                if 0 <= (t_term - dt) <= delta_t:
                    # Termination case 1
                    # Правильный случай прекращения напыления по уровню
                    term_cond = True
                    term_cond_case[j] = 1

                    # delta_t = nonloc_alg.calc_delta_t(dt, term_flux_lvl, set_up_pars.q_TR[j]) \
                    #     + norm_3sigma_rnd(rng, mean=set_up_pars.delay_time, sigma=set_up_pars.delay_time_sigma)
                    delta_t = t_term - dt

                    # шаг *изменения состояния природы*
                    cur_rate = set_up_pars.rates[j] + norm_3sigma_rnd(rng, sigma=set_up_pars.rates_sigmas[j])
                    des_act.increase_layer_thickness(j, delta_t * cur_rate)

                    time_list[j].append(time_list[j][-1] + delta_t)
                    d_j_act_t[j].append(des_act.d[j])
                    flux_act[j].append(des_act.calc_flux(set_up_pars.waves[j], q_TR=set_up_pars.q_TR[j], layer=j,
                                                         backside=set_up_pars.backside, q_subs=q_subs, save_M=save_M))
                    flux_meas[j].append(flux_act[j][-1] + norm_3sigma_rnd(rng, sigma=set_up_pars.meas_sigmas[j]))

                elif dt > t_term:
                    # Termination case 2
                    # На очередном шаге было обнаружено, что произошло перепыление
                    term_cond = True
                    term_cond_case[j] = 2

                elif des_act.d[j] > 1.25 * des_th.d[j]:
                    # Termination case 3
                    # Экстренное прекращение напыления --- остальные критерии не сработали
                    term_cond = True
                    term_cond_case[j] = 3

                if term_cond:
                    if set_up_pars.delay_time_sigma != 0.0:
                        # Ошибка времени закрытия заслонки
                        shutter_delay_dt = norm_3sigma_rnd(rng, mean=set_up_pars.delay_time, sigma=set_up_pars.delay_time_sigma)
                        # На этом шаге прирост толщины может быть и положительным и отрицательным
                        des_act.increase_layer_thickness(j, shutter_delay_dt * cur_rate)
                        d_j_act_t[j][-1] = des_act.d[j]
                    break
        else:
            raise NameError(f'Overflow simulation error on seed = {rnd_seed}')

    wavelength = [set_up_pars.waves[j].wavelength for j in range(1, des_th.N + 1)]
    return SimInfo(des_th.d, des_act.d, time_list, flux_meas, term_cond_case, wavelength,
                   rnd_seed=rnd_seed, d_j_act_t=d_j_act_t, set_up_pars=set_up_pars, des=des_th)


# def measurement_simulation(sim_info: SimInfo, waves):
#     rnd_seed = sim_info.rnd_seed
#     rng = np.random.default_rng(rnd_seed)
#
#     des_th = sim_info.des
#     set_up_pars = sim_info.set_up_pars
#     # Для производительности определим и выделим заранее достаточное кол-во памяти массивам
#     len_sim_list = num_step_layer_estimation(des_th, set_up_pars)
#     max_steps = len_sim_list[1]
#
#     time_list = [[0.0] for _ in range(des_th.N + 1)]
#     d_j_act_t = [[0.0] for _ in range(des_th.N + 1)]
#     flux_meas = [[0.0] for _ in range(des_th.N + 1)]
#     flux_act = [[0.0] for _ in range(des_th.N + 1)]
#
#     # Информация о дизайне в текущий для симуляции момент времени
#     des_act = copy.deepcopy(des_th)
#     # des_act.d = np.zeros(des_th.N + 1, dtype=float)
#     des_act.d = (des_th.N + 1) * [0.0]
#     # term_cond_case = np.zeros(des_th.N + 1, dtype=float)
#     term_cond_case = (des_th.N + 1) * [0]
#
#     for j in range(1, des_th.N + 1):
#         term_cond = False
#
#         # в конце напыления шаг по времени может быть меньше tau
#         # поэтому введём переменную delta_t означающую размер интервала времени для выполняемого шага
#         delta_t = set_up_pars.tau
#         # dt --- времени прошло с начала напыления слоя
#         dt = 0.0
#         expected_interval = False
#
#         # шаг *измерения и анализа*
#         # lsc = 0  # layer_step_counter
#         time_list[j][0] = time_list[j - 1][-1]
#         flux_act[j][0] = des_act.calc_flux(set_up_pars.waves[j], q_TR=set_up_pars.q_TR[j], layer=j,
#                                            backside=set_up_pars.backside)
#         flux_meas[j][0] = flux_act[j][-1] + norm_3sigma_rnd(rng, sigma=set_up_pars.meas_sigmas[j])
#
#         for lsc in range(1, max_steps):
#
#             # шаг *изменения состояния природы*
#             dt += delta_t
#             cur_rate = set_up_pars.rates[j] + norm_3sigma_rnd(rng, sigma=set_up_pars.rates_sigmas[j])
#             des_act.increase_layer_thickness(j, delta_t * cur_rate)
#
#             # шаг *измерения и анализа*
#             time_list[j].append(time_list[j][-1] + delta_t)
#             d_j_act_t[j].append(des_act.d[j])
#             flux_act[j].append(des_act.calc_flux(set_up_pars.waves[j], q_TR=set_up_pars.q_TR[j], layer=j,
#                                                  backside=set_up_pars.backside))
#             flux_meas[j].append(flux_act[j][-1] + norm_3sigma_rnd(rng, sigma=set_up_pars.meas_sigmas[j]))
#
#                 t_term = nonloc_alg.calc_t(d_th=des_th.d[j], lvl=term_flux_lvl, q_TR=set_up_pars.q_TR[j])
#
#                 if 0 <= (t_term - dt) <= delta_t:
#                     # Termination case 1
#                     # Правильный случай прекращения напыления по уровню
#                     term_cond = True
#                     term_cond_case[j] = 1
#
#                     # delta_t = nonloc_alg.calc_delta_t(dt, term_flux_lvl, set_up_pars.q_TR[j]) \
#                     #     + norm_3sigma_rnd(rng, mean=set_up_pars.delay_time, sigma=set_up_pars.delay_time_sigma)
#                     delta_t = t_term - dt
#
#                     # шаг *изменения состояния природы*
#                     cur_rate = set_up_pars.rates[j] + norm_3sigma_rnd(rng, sigma=set_up_pars.rates_sigmas[j])
#                     des_act.increase_layer_thickness(j, delta_t * cur_rate)
#
#                     time_list[j].append(time_list[j][-1] + delta_t)
#                     d_j_act_t[j].append(des_act.d[j])
#                     flux_act[j].append(des_act.calc_flux(set_up_pars.waves[j], q_TR=set_up_pars.q_TR[j], layer=j,
#                                                          backside=set_up_pars.backside))
#                     flux_meas[j].append(flux_act[j][-1] + norm_3sigma_rnd(rng, sigma=set_up_pars.meas_sigmas[j]))
#
#                 elif dt > t_term:
#                     # Termination case 2
#                     # На очередном шаге было обнаружено, что произошло перепыление
#                     term_cond = True
#                     term_cond_case[j] = 2
#
#                 elif des_act.d[j] > 1.25 * des_th.d[j]:
#                     # Termination case 3
#                     # Экстренное прекращение напыления --- остальные критерии не сработали
#                     term_cond = True
#                     term_cond_case[j] = 3
#
#                 if term_cond:
#                     if set_up_pars.delay_time_sigma != 0.0:
#                         # Ошибка времени закрытия заслонки
#                         shutter_delay_dt = norm_3sigma_rnd(rng, mean=set_up_pars.delay_time, sigma=set_up_pars.delay_time_sigma)
#                         # На этом шаге прирост толщины может быть и положительным и отрицательным
#                         des_act.increase_layer_thickness(j, shutter_delay_dt * cur_rate)
#                         d_j_act_t[j][-1] = des_act.d[j]
#                     break
#         else:
#             raise NameError(f'Overflow simulation error on seed = {rnd_seed}')
#
#     wavelength = [set_up_pars.waves[j].wavelength for j in range(1, des_th.N + 1)]
#     return SimInfo(des_th.d, des_act.d, time_list, flux_meas, term_cond_case, wavelength, d_j_act_t=d_j_act_t)








