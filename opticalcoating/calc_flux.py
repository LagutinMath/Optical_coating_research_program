from functools import singledispatchmethod, total_ordering
import numpy as np
import scipy

@total_ordering
class Wave:
    __slots__ = ('wavelength', 'polarisation', 'angle')

    @singledispatchmethod
    def __init__(self, wavelength, polarisation='S', angle=0.):
        """Incident wave class
        :param wavelength: wavelength in nm
        :param polarisation: 'S' or 'P' polarisation
        :param angle: angle in radians
        """
        self.wavelength = wavelength
        self.polarisation = polarisation
        self.angle = angle

    @__init__.register(tuple)
    def _from_tuple(self, wave_tuple):
        self.wavelength = wave_tuple[0]
        self.polarisation = wave_tuple[1] if len(wave_tuple) > 1 and wave_tuple[1] is not None else 'S'
        self.angle = wave_tuple[2] if len(wave_tuple) > 2 and wave_tuple[2] is not None else 0

    @classmethod
    def from_tuple(cls, wave_tuple: tuple):
        wavelength = wave_tuple[0]
        polarisation = wave_tuple[1] if len(wave_tuple) > 1 and wave_tuple[1] is not None else 'S'
        angle = wave_tuple[2] if len(wave_tuple) > 2 and wave_tuple[2] is not None else 0
        return cls(wavelength, polarisation, angle)

    def __float__(self):
        return float(self.wavelength)

    def __repr__(self):
        return f"Wave({self.wavelength}, '{self.polarisation}', {self.angle})"

    def __eq__(self, other):
        if isinstance(other, Wave):
            return (self.wavelength == other.wavelength and
                    self.polarisation == other.polarisation and
                    self.angle == other.angle)
        return NotImplemented

    def __lt__(self, other):
        if isinstance(other, Wave):
            if self.angle > other.angle: return True
            if self.polarisation != other.polarisation and other.polarisation == 'P': return True
            return self.wavelength < other.wavelength
        return NotImplemented


def gauss(flux_list):
    """Численное пятиточечное интегрирование по формулам Гаусса-Кристоффеля с гауссовой весовой функцией"""
    return (0.011257411327720682 * flux_list[0] + 0.22207592200561266 * flux_list[1] +
            0.5333333333333333 * flux_list[2] +
            0.22207592200561266 * flux_list[3] + 0.011257411327720682 * flux_list[4])


def sigma(FHWM):
    # FHWM = 2 Sqrt[2 Log[2]] * sigma
    return 0.5887050112577373 * FHWM


def sp_mat_mult(A, B):  # special matrix multiplication
    m00 = A[0][0] * B[0][0] - A[0][1] * B[1][0]
    m01 = A[0][0] * B[0][1] + A[0][1] * B[1][1]
    m10 = A[1][0] * B[0][0] + A[1][1] * B[1][0]
    m11 = A[1][1] * B[1][1] - A[1][0] * B[0][1]
    return [[m00, m01], [m10, m11]]


def calc_flux(des, wv, *, q_subs=True, backside=False, q_percent=False, n_a=1, q_TR='R', layer=None, save_M=False, width=None, width_form='rect'):
    if layer is not None:
        w_num = des.witness_num(layer)
    else:
        w_num = des.w_cur

    # упрощенная процедура вычисления с использованием вычисленной ранее хар-ой матрицы
    if save_M:
        if des.j_changed:
            des.cur_M = des.temp_M
            des.j_changed = False
        M = des.cur_M
        j = des.cur_j
        if wv.angle == 0:
            q_j = n_j = des.n(j, wv)
            phi_j = 2.0 * np.pi * n_j * des.d[j] / wv.wavelength
        else:
            n_j = des.n(j, wv)
            gamma_j = np.arcsin(np.sin(wv.angle) * n_a / n_j)
            q_j = n_j * np.cos(gamma_j) if wv.polarisation == 'S' else n_j / np.cos(gamma_j)
            phi_j = 2.0 * np.pi * n_j * des.d[j] * np.cos(gamma_j) / wv.wavelength

        cos_phi = np.cos(phi_j)
        sin_phi = np.sin(phi_j)

        M_j_00 = cos_phi
        M_j_01 = sin_phi / q_j
        M_j_10 = sin_phi * q_j
        M_j_11 = cos_phi

        M = sp_mat_mult([[M_j_00, M_j_01], [M_j_10, M_j_11]], M)
        des.temp_M = M
    else:
        M = [[1.0, 0.0], [0.0, 1.0]]
        for j in des.witnesses[w_num]:
            if wv.angle == 0:
                q_j = n_j = des.n(j, wv)
                phi_j = 2.0 * np.pi * n_j * des.d[j] / wv.wavelength
            else:
                n_j = des.n(j, wv)
                gamma_j = np.arcsin(np.sin(wv.angle) * n_a / n_j)
                q_j = n_j * np.cos(gamma_j) if wv.polarisation == 'S' else n_j / np.cos(gamma_j)
                phi_j = 2.0 * np.pi * n_j * des.d[j] * np.cos(gamma_j) / wv.wavelength

            cos_phi = np.cos(phi_j)
            sin_phi = np.sin(phi_j)

            M_j_00 = cos_phi
            M_j_01 = sin_phi / q_j
            M_j_10 = sin_phi * q_j
            M_j_11 = cos_phi

            M = sp_mat_mult([[M_j_00, M_j_01], [M_j_10, M_j_11]], M)

    n_s = des.n(0, wv)
    if wv.angle == 0:
        q_s = n_s
        q_a = n_a
    else:
        gamma_s = np.arcsin(np.sin(wv.angle) * n_a / n_s)
        q_s = n_s * np.cos(gamma_s) if wv.polarisation == 'S' else n_s / np.cos(gamma_s)
        q_a = n_a * np.cos(wv.angle) if wv.polarisation == 'S' else n_a / np.cos(wv.angle)

    T_1 = 4.0 * q_a * q_s / ((M[0][0] * q_a + M[1][1] * q_s) ** 2 + (M[1][0] + M[0][1] * q_a * q_s) ** 2)
    R_1 = (((M[0][0] * q_a - M[1][1] * q_s) ** 2 + (M[0][1] * q_a * q_s - M[1][0]) ** 2) /
           ((M[0][0] * q_a + M[1][1] * q_s) ** 2 + (M[0][1] * q_a * q_s + M[1][0]) ** 2))
    if q_subs:
        T_2 = 4.0 * q_a * q_s / (q_a + q_s) ** 2
        R_2 = ((q_a - q_s) / (q_a + q_s)) ** 2
        if backside:
            T_1, T_2 = T_2, T_1
            R_1, R_2 = R_2, R_1
        # формула для коэффициента поглощения верна для нормального падения
        sgm = np.exp(-4 * np.pi * des.xi(0, wv) * des.d[0] / wv.wavelength)
        T = sgm * T_1 * T_2 / (1.0 - R_1 * R_2 * sgm**2)
        R = (R_1 + (sgm**2 * R_2 * T_1 * T_1) / (1.0 - R_1 * R_2 * sgm**2))
    else:
        T, R = T_1, R_1

    if width is not None:
        if width_form in 'rect':
            flux_in_wv = lambda x: calc_flux(des, Wave(x, wv.polarisation, wv.angle), q_subs=q_subs,
                                             backside=backside, q_percent=q_percent, n_a=n_a,
                                             q_TR=q_TR, layer=layer, save_M=save_M)
            a, b = wv.wavelength - 0.5 * width, wv.wavelength + 0.5 * width
            return scipy.integrate.romberg(flux_in_wv, a, b, tol=0.0001) / width
        elif width_form in 'gauss':
            x_list = (-2.8569700138728056, -1.355626179974266, 0.0,
                   1.355626179974266, 2.8569700138728056)
            wvs = [Wave(wv.wavelength + sigma(width) * x) for x in x_list]
            flux_list = [calc_flux(des, wave, q_subs=q_subs, backside=backside, q_percent=q_percent, n_a=n_a,
                                   q_TR=q_TR, layer=layer, save_M=save_M) for wave in wvs]
            return gauss(flux_list)

    if q_TR == 'both':
        return (100.0 * T, 100.0 * R) if q_percent else (T, R)
    elif q_TR == 'T':
        return 100.0 * T if q_percent else T
    elif q_TR == 'R':
        return 100.0 * R if q_percent else R
    else:
        raise NameError('Wrong q_TR')


def grad_flux(des, wv, *, n_a=1, q_TR='R'):
    M = [[1.0, 0.0], [0.0, 1.0]]
    M_j = np.empty((des.N + 1, 2, 2), dtype=complex)
    D_j = np.empty((des.N + 1, 2, 2), dtype=complex)
    for j in range(1, des.N + 1):
        if wv.angle == 0:
            q_j = n_j = des.n(j, wv)
            phi_j = 2.0 * np.pi * n_j * des.d[j] / wv.wavelength
            cos_gamma_j = 1.0
        else:
            n_j = des.n(j, wv)
            gamma_j = np.arcsin(np.sin(wv.angle) * n_a / n_j)
            cos_gamma_j = np.cos(gamma_j)
            q_j = n_j * cos_gamma_j if wv.polarisation == 'S' else n_j / cos_gamma_j
            phi_j = 2.0 * np.pi * n_j * des.d[j] * cos_gamma_j / wv.wavelength

        cos_phi = np.cos(phi_j)
        sin_phi = np.sin(phi_j)

        M_j_00 = cos_phi
        M_j_01 = sin_phi / q_j
        M_j_10 = sin_phi * q_j
        M_j_11 = cos_phi

        M_j[j] = np.array([[M_j_00, 1j * M_j_01], [1j * M_j_10, M_j_11]])
        D_j[j] = (2.0 * np.pi * cos_gamma_j * n_j / wv.wavelength) * np.array(
            [[- sin_phi, (1j / q_j) * cos_phi], [1j * q_j * cos_phi, - sin_phi]])

        M = sp_mat_mult([[M_j_00, M_j_01], [M_j_10, M_j_11]], M)

    n_s = des.n(0, wv)
    if wv.angle == 0:
        q_s = n_s
        q_a = n_a
    else:
        gamma_s = np.arcsin(np.sin(wv.angle) * n_a / n_s)
        q_s = n_s * np.cos(gamma_s) if wv.polarisation == 'S' else n_s / np.cos(gamma_s)
        q_a = n_a * np.cos(wv.angle) if wv.polarisation == 'S' else n_a / np.cos(wv.angle)

    t = 2 * q_a / (q_a * M[0][0] + q_s * M[1][1] + 1j * (q_s * q_a * M[0][1] + M[1][0]))
    conj_t = np.conjugate(t)
    Psi_t = (- t ** 2 / (2 * q_a)) * np.array([[q_a, 1.0], [q_s * q_a, q_s]])

    grad_T = np.empty(des.N + 1, dtype=float)
    for j in range(1, des.N + 1):
        M_prod = prod(M_j, des.N, j + 1)
        Psi_prod = prod(M_j, j - 1, 1) @ Psi_t
        dt_j = np.trace(M_prod @ D_j[j] @ Psi_prod)
        grad_T[j] = 2 * q_s / q_a * np.real(conj_t * dt_j)

    if q_TR == 'T':
        return grad_T
    else:
        # case q_TR == 'R'
        return - grad_T


def prod(A, start, end):
    res = np.eye(2)

    if end > len(A) - 1:
        return res
    elif start < 1:
        return res

    for j in range(start, end - 1, -1):
        res = res@A[j]
    return res