# import os.path
# import numpy as np
# import math
from design_class import *
import copy
from os import listdir
from matplotlib.animation import FuncAnimation


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
                 rates_sigmas, rates_syst, delay_time, delay_time_sigma, r_index_syst, r_index_sigmas):
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
        self.delay_time_sigma = delay_time_sigma
        self.r_index_syst = r_index_syst
        self.r_index_sigmas = r_index_sigmas
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


def norm_3sigma_rnd(rng, *, mean=0.0, sigma=0.0):
    """
    Функция возвращает нормально распределенную случайную величину,
    с ненулевым мат.ожиданием, которая модифицирована так, чтобы
    значение не выходило за границы интервала [-3*sigma, 3*sigma]
    """
    rnd_val = rng.normal(mean, sigma)
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

    def q_expected_interval(self, dt, d_th, *, theta_th):
        """
        Сообщает находится ли в момент времени dt
        процесс напыления в промежутке между теми экстремумами,
        где должен прекратиться процесс напыления
        :param dt: время прошедшее с начала напыления слоя
        :param d_th: толщина теоретического дизайна
        :param theta_th: оценка теоретического theta на слое
        """
        phi_act = self.n * self.r * dt * 2 * np.pi / self.wavelength
        phi_th = self.n * d_th * 2 * np.pi / self.wavelength

        if num_cutoffs(theta_th, phi_act) == num_cutoffs(theta_th, phi_th):
            return True
        else:
            return False

    def prev_extr(self, dt, q_TR):
        phi_act = self.n * self.r * dt * 2 * np.pi / self.wavelength
        return self.coef.prev_extr(phi_act, q_TR)

    def q_pass_term_flux_lvl(self, dt, term_flux_lvl, q_TR):
        flux_cur = self.flux(self.r * dt, q_TR)
        flux_prev_extr = self.prev_extr(dt, q_TR)
        return ((flux_prev_extr - term_flux_lvl) * (flux_cur - term_flux_lvl)) < 0.0

    def calc_delta_t(self, dt, term_flux_lvl, q_TR, tau=2.0):
        phi_act = self.n * self.r * dt * 2 * np.pi / self.wavelength
        if q_TR == 'R':
            temp = (1.0 / (term_flux_lvl - 1.0) - self.coef.gamma()) / self.coef.D()
        else:
            temp = (1.0 / (- term_flux_lvl) - self.coef.gamma()) / self.coef.D()
        # + solution
        phi_lvl_p = 0.5 * (math.acos(temp) - self.coef.theta())
        # - solution
        phi_lvl_n = - 0.5 * (math.acos(temp) + self.coef.theta())
        # if q_TR == 'R':
        #     if (term_flux_lvl - self.prev_extr(dt, q_TR) >= 0.0) and (math.sin(2 * phi_lvl_p + self.coef.theta()) >= 0.0):
        #         q_p_sol = True
        #     elif (term_flux_lvl - self.prev_extr(dt, q_TR) < 0.0) and (math.sin(2 * phi_lvl_p + self.coef.theta()) < 0.0):
        #         q_p_sol = True
        #     else:
        #         q_p_sol = False
        # else:
        #     if (term_flux_lvl - self.prev_extr(dt, q_TR) >= 0.0) and (math.sin(2 * phi_lvl_p + self.coef.theta()) < 0.0):
        #         q_p_sol = True
        #     elif (term_flux_lvl - self.prev_extr(dt, q_TR) < 0.0) and (math.sin(2 * phi_lvl_p + self.coef.theta()) >= 0.0):
        #         q_p_sol = True
        #     else:
        #         q_p_sol = False
        #
        # phi_lvl = phi_lvl_p if q_p_sol else phi_lvl_n
        #
        # delta_phi = (phi_lvl - phi_act) - 0.5 * np.pi * math.floor(2 * (phi_lvl - phi_act)/np.pi)
        # delta_t = delta_phi * self.wavelength / (self.n * self.r * 2 * np.pi)

        # Грубое решение о выборе корня
        phi_lvl = phi_lvl_p
        delta_phi = (phi_lvl - phi_act) - np.pi * math.floor((phi_lvl - phi_act)/np.pi)
        delta_t_p = delta_phi * self.wavelength / (self.n * self.r * 2 * np.pi)

        phi_lvl = phi_lvl_n
        delta_phi = (phi_lvl - phi_act) - np.pi * math.floor((phi_lvl - phi_act) / np.pi)
        delta_t_n = delta_phi * self.wavelength / (self.n * self.r * 2 * np.pi)

        if (delta_t_p >= 0.) and (delta_t_p <= tau):
            return delta_t_p
        elif (delta_t_n >= 0.) and (delta_t_n <= tau):
            return delta_t_n
        else:
            return 0.0


def find_file_name(obj_name: str, ext='.json'):
    """
    Ищет и возвращает свободное имя в папке *имя объекта*s
    Пример: obj_name = 'Simulation', первое возвращенное имя будет:
    'Simulations/Simulation001.json'
    """
    dir_file_names = listdir(obj_name + 's/')
    indx = 1
    file_name = obj_name + 's/' + obj_name + str(indx).zfill(3) + ext

    while indx < 1000:
        file_name = obj_name + 's/' + obj_name + str(indx).zfill(3) + ext
        is_name_in_dir = False
        for name in dir_file_names:
            is_name_in_dir = ((obj_name + 's/' + name) == file_name)
            if is_name_in_dir:
                break
        if is_name_in_dir:
            indx += 1
            continue
        else:
            break

    return file_name


class SimInfo:
    def __init__(self, des_th_d, des_act_d, time_list_res, flux_meas_res, term_cond_case):
        self.N_layers = len(des_th_d) - 1
        self.d_th = des_th_d
        self.d_act = des_act_d
        self.time_list = time_list_res
        self.flux_meas = flux_meas_res
        self.term_cond_case = term_cond_case
        # self.errors_d = des_act_d - des_th_d
        self.errors_d = [d_act - d_th for (d_act, d_th) in list(zip(*[des_act_d, des_th_d]))]

    def make_dict(self):
        len_list = len(self.time_list)
        time_list = len_list * [[]]
        flux_meas = len_list * [[]]
        for j in range(1, len_list):
            time_list[j] = self.time_list[j].tolist()
            flux_meas[j] = self.flux_meas[j].tolist()
        sim_dict = {'time_list': time_list, 'flux_meas': flux_meas,
                    'term_cond_case': self.term_cond_case, 'errors_d': self.errors_d}
        return sim_dict

    def save(self):
        file_name = find_file_name('Simulation')

        with open(file_name, 'w') as file:
            json.dump(self.make_dict(), file, indent=3)
            file.close()

    def animation(self, j=1):
        """
        Анимация напыления j-ого слоя
        """
        x_len = len(self.time_list[j])
        x_min = self.time_list[j][0]
        x_max = self.time_list[j][x_len - 1]

        fig, ax = plt.subplots()
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(0., 1.)
        line, = ax.plot(0., 0.)

        x_data = []
        y_data = []

        def animation_frame(i: int):
            """
            Отрисовка линии до i_end точки на j-ом слое
            """
            x_data.append(self.time_list[j][i])
            y_data.append(self.flux_meas[j][i])

            line.set_xdata(x_data)
            line.set_ydata(y_data)
            return line,

        animation = FuncAnimation(fig, func=animation_frame, frames=range(x_len), interval=1)
        plt.show()


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

        self.nonloc_coef = (des_th.N + 1) * [[]]

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
            self.nonloc_coef[j_cur] = [abg.alpha(), abg.beta(), abg.gamma(), abg.D(), abg.theta()]

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
    max_steps = len_sim_list[1]

    # time_list = np.empty(len_sim_list, dtype=float)
    # flux_meas = np.empty(len_sim_list, dtype=float)
    # flux_act = np.empty(len_sim_list, dtype=float)

    time_list = [[0.0] for _ in range(des_th.N + 1)]
    flux_meas = [[0.0] for _ in range(des_th.N + 1)]
    flux_act = [[0.0] for _ in range(des_th.N + 1)]

    # Информация о дизайне в текущий для симуляции момент времени
    des_act = copy.deepcopy(des_th)
    # des_act.d = np.zeros(des_th.N + 1, dtype=float)
    des_act.d = (des_th.N + 1) * [0.0]

    # term_cond_case = np.zeros(des_th.N + 1, dtype=float)
    term_cond_case = (des_th.N + 1) * [0]

    for j in range(1, des_th.N + 1):
        des_act.n_fix[0][j] = des_act.n(j, set_up_pars.waves[j])\
                * (1 + norm_3sigma_rnd(rng, mean=set_up_pars.r_index_syst[j], sigma=set_up_pars.r_index_sigmas[j]))
        # На новом слое получем теоретическое theta (можно просто оценить, но мы получим точно)
        nonloc_coef_th = MonochromStrategyInfo(des_act, set_up_pars).nonloc_coef[j]
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
        flux_act[j][0] = calc_flux(des_act, set_up_pars.waves[j], q_TR=set_up_pars.q_TR[j])
        flux_meas[j][0] = flux_act[j][-1] + norm_3sigma_rnd(rng, sigma=set_up_pars.meas_sigmas[j])
        nonloc_alg.refresh(dt, flux_meas[j][-1], set_up_pars.q_TR[j])

        # напыление слоя
        # while not term_cond:
        for lsc in range(1, max_steps):
            # шаг *изменения состояния природы*
            dt += delta_t
            cur_rate = set_up_pars.rates[j] + norm_3sigma_rnd(rng, sigma=set_up_pars.rates_sigmas[j])
            des_act.increase_layer_thickness(j, delta_t * cur_rate)

            # шаг *измерения и анализа*
            time_list[j].append(time_list[j][-1] + delta_t)
            flux_act[j].append(calc_flux(des_act, set_up_pars.waves[j], q_TR=set_up_pars.q_TR[j]))
            flux_meas[j].append(flux_act[j][-1] + norm_3sigma_rnd(rng, sigma=set_up_pars.meas_sigmas[j]))
            nonloc_alg.refresh(dt, flux_meas[j][-1], set_up_pars.q_TR[j])

            if not nonloc_alg.exist_coef:
                continue
            elif not expected_interval:
                if des_act.d[j] < 0.7 * des_th.d[j]:
                    continue
                else:
                    expected_interval = nonloc_alg.q_expected_interval(dt, des_th.d[j], theta_th=theta_th)
            elif expected_interval:
                if term_algs[j] == 'Elimination':
                    term_flux_lvl = nonloc_alg.flux(des_th.d[j], set_up_pars.q_TR[j])
                elif term_algs[j] == 'Quasiswing':
                    turn = nonloc_alg.prev_extr(dt, set_up_pars.q_TR[j])
                    ampl = nonloc_alg.coef.ampl()
                    term_flux_lvl = turn - (str_info.prev_extr[j] - str_info.term[j]) * (ampl / str_info.ampl[j])
                elif term_algs[j] == 'Pseudoswing':
                    # Не готово
                    term_flux_lvl = 0.0
                elif term_algs[j] == 'None':
                    term_flux_lvl = str_info.term[j]
                else:
                    print('Uncorrect termination algoritm for layer', j)
                    break  # место в будущем для практики выкидывания ошибок

                if nonloc_alg.q_pass_term_flux_lvl(dt + delta_t, term_flux_lvl, set_up_pars.q_TR[j]):
                    # Termination case 1
                    # Правильный случай прекращения напыления по уровню
                    term_cond = True
                    term_cond_case[j] = 1

                    delta_t = nonloc_alg.calc_delta_t(dt, term_flux_lvl, set_up_pars.q_TR[j]) \
                        + norm_3sigma_rnd(rng, mean=set_up_pars.delay_time, sigma=set_up_pars.delay_time_sigma)


                    # шаг *изменения состояния природы*
                    cur_rate = set_up_pars.rates[j] + norm_3sigma_rnd(rng, sigma=set_up_pars.rates_sigmas[j])
                    des_act.increase_layer_thickness(j, delta_t * cur_rate)

                    if j < des_th.N:
                        time_list[j].append(time_list[j][-1] + delta_t)
                        flux_act[j].append(calc_flux(des_act, set_up_pars.waves[j], q_TR=set_up_pars.q_TR[j]))
                        flux_meas[j].append(flux_act[j][-1] + norm_3sigma_rnd(rng, sigma=set_up_pars.meas_sigmas[j]))

                elif nonloc_alg.q_pass_term_flux_lvl(dt, term_flux_lvl, set_up_pars.q_TR[j]):
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
                    break
            else:
                print('Impossible! (expected_interval block)')
        else:
            print('Overflow simulation error on seed =', rnd_seed)

    return SimInfo(des_th.d, des_act.d, time_list, flux_meas, term_cond_case)


class StatInfo:
    def __init__(self, des_th: Design, term_algs, set_up_pars, err_list, start_rnd_seed):
        self.err_list = err_list
        self.start_rnd_seed = start_rnd_seed
        self.des = des_th
        self.set_up_pars = set_up_pars
        self.term_algs = term_algs

        self.N_sim = len(err_list)

    def make_dict(self):
        stat_dict = {'design': self.des.name, 'start_rnd_seed': self.start_rnd_seed, 'error list': self.err_list,
                     'term_algs': self.term_algs}
        return stat_dict

    def save(self):
        file_name = find_file_name('Statistic')

        with open(file_name, 'w') as file:
            json.dump(self.make_dict(), file, indent=3)
            file.close()

    def save_plain_txt(self):
        file_name = find_file_name('Statistic', ext='.txt')

        n = len(self.err_list)
        m = len(self.err_list[0])

        with open(file_name, 'w') as file:
            for i in range(n):
                for j in range(m):
                    print(self.err_list[i][j], end='\t', file=file)
                print('', file=file)
            file.close()


def collect_statistics(des_th, term_algs, set_up_pars, N_sim, start_rnd_seed=10000000):

    err_list = N_sim * [[]]

    for x in range(N_sim):
        rnd_seed = start_rnd_seed + x
        res = simulation(des_th, term_algs, set_up_pars, rnd_seed)
        err_list[x] = res.errors_d[1:des_th.N + 1]

    return StatInfo(des_th, term_algs, set_up_pars, err_list, start_rnd_seed)



