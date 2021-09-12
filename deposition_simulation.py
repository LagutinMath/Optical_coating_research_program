# import os.path
# import numpy as np
# import math
from design_class import *
import copy
from os import listdir


class SetUpParameters:
    """
    Класс содержащий в себе описание
    1. Параметров установки заданных пользователем, таких как:
        - длина волны мониторинга,
        - наблюдается коэф.отражения или пропускания
        - промежуток времени через который проводятся измерения
    2. Информацию о несовершенстве установки
    """

    def __init__(self, waves, q_TR, tau, rates, mono_width, meas_sigmas, meas_syst,
                 rates_sigmas, rates_syst, delay_time):
        self.waves = waves
        self.q_TR = q_TR
        self.tau = tau
        self.rates = rates
        self.mono_width = mono_width
        self.meas_sigmas = meas_sigmas
        self.meas_syst = meas_syst
        self.rates_sigmas = rates_sigmas
        self.rates_syst = rates_syst
        self.delay_time = delay_time
        # N --- кол-во слоёв
        # в массивах индекс соответствует номеру слоя
        # индекс 0 незадействован для rates и пр.
        self.N = len(waves) - 1

    def strategy_spreadsheet_plot(self):
        """
        Построение кривой монохроматического мониторинга
        """
        pass


def num_step_estimation(des: Design, set_up_pars: SetUpParameters):
    """
    Вычисление оценки кол-ва шагов в предстоящей симуляции
    """
    r = set_up_pars.rates
    tau = set_up_pars.tau
    total_time = 0.0
    for j in range(1, des.N + 1):
        total_time += des.d[j] / r[j]
    return 10 * math.ceil(total_time / tau)


def num_step_layer_estimation(des: Design, set_up_pars: SetUpParameters):
    """
    Возвращает лист нужный для генерации пустых списков для предстоящей симуляции
    первое число кол-во слоёв (кол-во строк)
    второе число кол-во шагов на одном слое (кол-во столбцов)
    кол-во шагов возвращает с 10-ти кратным запасом
    """
    r = set_up_pars.rates
    tau = set_up_pars.tau
    max_step = 0
    for j in range(1, des.N + 1):
        max_step = max(max_step, des.d[j] / r[j])
    max_step = math.ceil(max_step / tau)
    return des.N + 1, 10 * max_step


def norm_3sigma_rnd(rng, sigma):
    """
    Функция возвращает нормально распределенную случайную величину,
    с нулевым мат.ожиданием, которая модифицирована так, чтобы
    значение не выходило за границы интервала [-3*sigma, 3*sigma]
    """
    rnd_val = rng.normal(0.0, sigma)
    if rnd_val > 3.0 * sigma:
        rnd_val = 3.0 * sigma
    elif rnd_val < -3.0 * sigma:
        rnd_val = -3.0 * sigma
    return rnd_val


# @njit(fastmath=True)
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


# @njit(fastmath=True)
def num_cutoffs(theta: float, phi: float):
    """
    Количество точек pi * n принадлежащих (theta, theta + 2 * phi), где n --- целое
    """
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
        """
        Нелокальные коэф. можно определить
        либо как abg = (alpha, beta, gamma),
        либо как Dthg = (D, theta, gamma)
        """
        if abg is None:
            self.abg = [0.0, 0.0, 0.0]
        else:
            self.abg = abg
        self.Dthg = Dthg

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
        return self.flux_max() - self.flux_min()

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


def theor_nonloccoef(A_re: float, A_im: float, n_j: float, n_s: float):
    """
    Вычисляет нелок.коэф-ты по адмиттансу (с учетом отражения от задней стороны подложки)
    :param A_re: Действительная часть адмиттанса
    :param A_im: Мнимая часть адмиттанса
    :param n_j: показатель преломления текущего слоя
    :param n_s: показатель преломления
    :return: Нелокальные коэффициенты
    """
    k_R = ((n_s - 1.) / (n_s + 1)) ** 2
    rho = math.sqrt(A_re ** 2 + A_im ** 2)

    alpha = (rho ** 2 - n_j ** 2) * (1. - n_j ** 2) / (8. * A_re * n_j ** 2)
    beta = (1 - n_j ** 2) * A_im / (4. * n_j * A_re)
    gamma = 0.5 + (1. / (k_R - 1.)) - ((rho ** 2 + n_j ** 2) * (1. + n_j ** 2) / (8. * A_re * n_j ** 2))

    return NonlocCoef([alpha, beta, gamma])


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
        """
        Возвращает значение flux, которое соответствует напыленной толщине thickness
        При подстановки d_j^{theor} это есть уровень остановки независимого метода контроля
        """
        d_phi = self.n * thickness * 2 * np.pi / self.wavelength
        if q_TR == 'R':
            ans = 1.0 + 1.0 / (self.coef.D() * math.cos(2 * d_phi + self.coef.theta()) + self.coef.gamma())
        else:
            ans = -1.0 / (self.coef.D() * math.cos(2 * d_phi + self.coef.theta()) + self.coef.gamma())
        return ans

    def q_expected_interval(self, dt, d_th):
        """
        Сообщает находится ли в момент времени dt
        процесс напыления в промежутке между теми экстремумами,
        где должен прекратиться процесс напыления
        :param dt: время прошедшее с начала напыления слоя
        :param d_th: толщина теоретического дизайна
        """
        phi_act = self.n * self.r * dt * 2 * np.pi / self.wavelength
        phi_th = self.n * d_th * 2 * np.pi / self.wavelength
        theta = self.coef.theta()
        if num_cutoffs(theta, phi_act) == num_cutoffs(theta, phi_th):
            return True
        else:
            return False

    # def q_pass_term_flux_lvl(self, dt, term_flux_lvl, q_TR):
    #     raw_phi_act = self.n * self.r * dt * 2 * np.pi / self.wavelength
    #     phi_act = raw_phi_act - 0.5 * np.pi * num_cutoffs(0.0, raw_phi_act)
    #     print('raw_phi_act =', raw_phi_act / np.pi)
    #     print('phi_act =', phi_act/np.pi)
    #     if q_TR == 'R':
    #         temp = (1.0 / (term_flux_lvl - 1.0) - self.coef.gamma()) / self.coef.D()
    #     else:
    #         temp = (1.0 / (- term_flux_lvl) - self.coef.gamma()) / self.coef.D()
    #     raw_phi_lvl = 0.5 * (math.acos(temp) - self.coef.theta())
    #     phi_lvl = raw_phi_lvl - 0.5 * np.pi * num_cutoffs(0.0, raw_phi_lvl)
    #     print('raw_phi_lvl =', raw_phi_lvl / np.pi)
    #     print('phi_lvl =', phi_lvl / np.pi)
    #     return phi_lvl < phi_act

    def prev_extr(self, dt, q_TR):
        phi_act = self.n * self.r * dt * 2 * np.pi / self.wavelength
        return self.coef.prev_extr(phi_act, q_TR)

    def q_pass_term_flux_lvl(self, dt, term_flux_lvl, q_TR):
        flux_cur = self.flux(self.r * dt, q_TR)
        flux_prev_extr = self.prev_extr(dt, q_TR)
        return ((flux_prev_extr - term_flux_lvl) * (flux_cur - term_flux_lvl)) < 0.0

    def calc_delta_t(self, dt, term_flux_lvl, q_TR):
        phi_act = self.n * self.r * dt * 2 * np.pi / self.wavelength
        if q_TR == 'R':
            temp = (1.0 / (term_flux_lvl - 1.0) - self.coef.gamma()) / self.coef.D()
        else:
            temp = (1.0 / (- term_flux_lvl) - self.coef.gamma()) / self.coef.D()
        # + solution
        phi_lvl_p = 0.5 * (math.acos(temp) - self.coef.theta())
        # - solution
        phi_lvl_n = - 0.5 * (math.acos(temp) + self.coef.theta())
        if q_TR == 'R':
            if (term_flux_lvl - self.prev_extr(dt, q_TR) >= 0.0) and (math.sin(2 * phi_lvl_p + self.coef.theta()) >= 0.0):
                q_p_sol = True
            elif (term_flux_lvl - self.prev_extr(dt, q_TR) < 0.0) and (math.sin(2 * phi_lvl_p + self.coef.theta()) < 0.0):
                q_p_sol = True
            else:
                q_p_sol = False
        else:
            if (term_flux_lvl - self.prev_extr(dt, q_TR) >= 0.0) and (math.sin(2 * phi_lvl_p + self.coef.theta()) < 0.0):
                q_p_sol = True
            elif (term_flux_lvl - self.prev_extr(dt, q_TR) < 0.0) and (math.sin(2 * phi_lvl_p + self.coef.theta()) >= 0.0):
                q_p_sol = True
            else:
                q_p_sol = False

        phi_lvl = phi_lvl_p if q_p_sol else phi_lvl_n

        delta_phi = (phi_lvl - phi_act) - 0.5 * np.pi * math.floor(2 * (phi_lvl - phi_act)/np.pi)
        delta_t = delta_phi * self.wavelength / (self.n * self.r * 2 * np.pi)
        return delta_t


class SimInfo:
    def __init__(self, des_th_d, des_act_d, time_list_res, flux_meas_res, term_cond_case):
        self.N_layers = len(des_th_d) - 1
        self.d_th = des_th_d
        self.d_act = des_act_d
        self.time_list = time_list_res
        self.flux_meas = flux_meas_res
        self.term_cond_case = term_cond_case
        self.errors_d = des_act_d - des_th_d

    def make_dict(self):
        len_list = len(self.time_list)
        time_list = len_list * [[]]
        flux_meas = len_list * [[]]
        for j in range(1, len_list):
            time_list[j] = self.time_list[j].tolist()
            flux_meas[j] = self.flux_meas[j].tolist()
        sim_dict = {'time_list': time_list, 'flux_meas': flux_meas,
                    'term_cond_case': self.term_cond_case.tolist(), 'errors_d': self.errors_d.tolist()}
        return sim_dict

    def save(self):
        dir_file_names = listdir('Simulations/')
        indx = 1
        while indx < 1000:
            file_name = 'Simulations/Simulation' + str(indx).zfill(3) + '.json'
            is_name_in_dir = False
            for name in dir_file_names:
                is_name_in_dir = (('Simulations/' + name) == file_name)
                if is_name_in_dir:
                    break
            if is_name_in_dir:
                indx += 1
                continue
            else:
                break

        with open(file_name, 'w') as file:
            json.dump(self.make_dict(), file, indent=3)


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
    def __init__(self, des_th: Design, set_up_pars: SetUpParameters):
        """
        Содержит информацию о значениях энергетических коэф. T/R
        в начале напыления слоя, в конце напыления слоя и экстремумы,
        которые могут достигаться при достаточной толщине слоя
        для заданной стратегии монохроматического мониторинга
        в случае теоретического дизайна с однородными слоями
        """
        self.start = (des_th.N + 1) * [0.0]
        self.max = (des_th.N + 1) * [0.0]
        self.min = (des_th.N + 1) * [0.0]
        self.term = (des_th.N + 1) * [0.0]
        self.ampl = (des_th.N + 1) * [0.0]
        self.prev_extr = (des_th.N + 1) * [0.0]

        for j_cur in range(1, des_th.N + 1):
            wavelength = set_up_pars.waves[j_cur].wavelength
            n_s = des_th.n(0, set_up_pars.waves[j_cur])
            A_re = n_s
            A_im = 0.
            for j in range(1, j_cur):
                d = des_th.d[j]
                n = des_th.n(j, set_up_pars.waves[j_cur])
                A_re, A_im = increase_admittance(A_re, A_im, wavelength, d, n)

            n_j_cur = des_th.n(j_cur, set_up_pars.waves[j_cur])

            abg = theor_nonloccoef(A_re, A_im, n_j_cur, n_s)
            abg.calc_Dthg()

            if set_up_pars.q_TR[j_cur] == 'T':
                self.start[j_cur] = - 1. / (abg.alpha() + abg.gamma())
                self.max[j_cur] = - 1. / (abg.gamma() + abg.D())
                self.min[j_cur] = - 1. / (abg.gamma() - abg.D())
                self.ampl[j_cur] = self.max[j_cur] - self.min[j_cur]
            else:
                self.start[j_cur] = 1. + 1. / (abg.alpha() + abg.gamma())
                self.max[j_cur] = 1. + 1. / (abg.gamma() - abg.D())
                self.min[j_cur] = 1. + 1. / (abg.gamma() + abg.D())
                self.ampl[j_cur] = self.max[j_cur] - self.min[j_cur]

            d_j_cur = des_th.d[j_cur]
            phi = (2. * np.pi / wavelength) * d_j_cur * n_j_cur

            if set_up_pars.q_TR[j_cur] == 'T':
                self.term[j_cur] = - 1. / (abg.alpha() * math.cos(2. * phi) + abg.beta() * math.sin(2. * phi) + abg.gamma())
            else:
                self.term[j_cur] = 1. + 1. / (abg.alpha() * math.cos(2. * phi) + abg.beta() * math.sin(2. * phi) + abg.gamma())

            self.prev_extr[j_cur] = abg.prev_extr(phi, set_up_pars.q_TR[j_cur])


def simulation(des_th, term_algs, set_up_pars, rnd_seed=10000000):
    str_info = MonochromStrategyInfo(des_th, set_up_pars)
    rng = np.random.default_rng(rnd_seed)
    # Для производительности определим и выделим заранее достаточное кол-во памяти массивам
    len_sim_list = num_step_layer_estimation(des_th, set_up_pars)

    time_list = np.empty(len_sim_list, dtype=float)
    flux_meas = np.empty(len_sim_list, dtype=float)
    flux_act = np.empty(len_sim_list, dtype=float)

    time_list_res = (des_th.N + 1) * [[]]
    flux_meas_res = (des_th.N + 1) * [[]]

    # Информация о дизайне в текущий для симуляции момент времени
    des_act = copy.deepcopy(des_th)
    des_act.d = np.zeros(des_th.N + 1, dtype=float)

    term_cond_case = np.zeros(des_th.N + 1, dtype=float)

    time_list[1][0] = 0.0

    for j in range(1, des_th.N + 1):
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
        lsc = 0  # layer_step_counter
        flux_act[j][lsc] = calc_flux(des_act, set_up_pars.waves[j], q_TR=set_up_pars.q_TR[j])
        flux_meas[j][lsc] = flux_act[j, lsc] + norm_3sigma_rnd(rng, set_up_pars.meas_sigmas[j])
        nonloc_alg.refresh(dt, flux_meas[j][lsc], set_up_pars.q_TR[j])

        # напыление слоя
        while not term_cond:
            # шаг *изменения состояния природы*
            dt += delta_t
            cur_rate = set_up_pars.rates[j] + norm_3sigma_rnd(rng, set_up_pars.rates_sigmas[j])
            des_act.increase_layer_thickness(j, delta_t * cur_rate)

            # шаг *измерения и анализа*
            lsc += 1
            time_list[j][lsc] = time_list[j][lsc - 1] + delta_t
            flux_act[j][lsc] = calc_flux(des_act, set_up_pars.waves[j], q_TR=set_up_pars.q_TR[j])
            flux_meas[j][lsc] = flux_act[j, lsc] + norm_3sigma_rnd(rng, set_up_pars.meas_sigmas[j])
            nonloc_alg.refresh(dt, flux_meas[j][lsc], set_up_pars.q_TR[j])

            if not nonloc_alg.exist_coef:
                continue
            elif not expected_interval:
                if des_act.d[j] < 0.7 * des_th.d[j]:
                    continue
                else:
                    expected_interval = nonloc_alg.q_expected_interval(dt, des_th.d[j])

            if expected_interval:
                if term_algs[j] == 'Elimination':
                    term_flux_lvl = nonloc_alg.flux(des_th.d[j], set_up_pars.q_TR[j])
                elif term_algs[j] == 'Quasiswing':
                    turn = nonloc_alg.prev_extr(dt, set_up_pars.q_TR[j])
                    ampl = nonloc_alg.coef.ampl()
                    term_flux_lvl = turn - (str_info.prev_extr[j] - str_info.term[j]) * (ampl / str_info.ampl[j])
                elif term_algs[j] == 'Pseudoswing':
                    # Не готово
                    term_flux_lvl = 0.0
                else:  # term_algs[j] == 'None':
                    # Не готово
                    term_flux_lvl = 0.0

                if nonloc_alg.q_pass_term_flux_lvl(dt + delta_t, term_flux_lvl, set_up_pars.q_TR[j]):
                    # Termination case 1
                    term_cond = True
                    term_cond_case[j] = 1

                    delta_t = nonloc_alg.calc_delta_t(dt, term_flux_lvl, set_up_pars.q_TR[j])
                    if j < des_th.N:
                        time_list[j + 1][0] = time_list[j][lsc] + delta_t

                    # шаг *изменения состояния природы*
                    cur_rate = set_up_pars.rates[j] + norm_3sigma_rnd(rng, set_up_pars.rates_sigmas[j])
                    des_act.increase_layer_thickness(j, delta_t * cur_rate)

                elif nonloc_alg.q_pass_term_flux_lvl(dt, term_flux_lvl, set_up_pars.q_TR[j]):
                    # Termination case 2
                    term_cond = True
                    term_cond_case[j] = 2

                elif des_act.d[j] > 1.25 * des_th.d[j]:
                    # Termination case 3
                    term_cond = True
                    term_cond_case[j] = 3

                if term_cond:
                    time_list_res[j] = time_list[j][0:lsc]
                    flux_meas_res[j] = flux_meas[j][0:lsc]
                    break

    res = SimInfo(des_th.d, des_act.d, time_list_res, flux_meas_res, term_cond_case)
    return res


def collect_statistics(des_th, term_algs, set_up_pars, N_sim, start_rnd_seed=10000000):

    # p1 = multiprocessing.Process(target=simulation, args=((des_th, term_algs, set_up_pars, start_rnd_seed),))
    err_list = N_sim * [[]]

    for x in range(N_sim):
        rnd_seed = start_rnd_seed + x
        res = simulation(des_th, term_algs, set_up_pars, rnd_seed)
        err_list[x] = res.errors_d[1:des_th.N + 1].tolist()
    return err_list



