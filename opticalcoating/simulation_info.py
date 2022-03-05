from opticalcoating.save_data import find_file_name
import json

class SimInfo:
    def __init__(self, des_th_d, des_act_d, time_list_res, flux_meas_res, term_cond_case, wavelength):
        self.N_layers = len(des_th_d) - 1
        self.d_th = des_th_d
        self.d_act = des_act_d
        self.time_list = time_list_res
        self.flux_meas = flux_meas_res
        self.term_cond_case = term_cond_case
        self.errors_d = [d_act - d_th for (d_act, d_th) in list(zip(*[des_act_d, des_th_d]))]
        self.wavelength = wavelength


    def make_dict(self):
        sim_dict = {'time_list': self.time_list, 'flux_meas': self.flux_meas,
                    'wavelength': self.wavelength,
                    'actual thicnesses': self.d_act, 'actual n': self.n_act}
        return sim_dict


    def save(self):
        file_name = find_file_name('Simulation')

        with open(file_name, 'w') as file:
            json.dump(self.make_dict(), file, indent=3)
            file.close()


    # def animation(self, j=1):
    #     """Анимация напыления j-ого слоя"""
    #     x_len = len(self.time_list[j])
    #     x_min = self.time_list[j][0]
    #     x_max = self.time_list[j][x_len - 1]
    #
    #     fig, ax = plt.subplots()
    #     ax.set_xlim(x_min, x_max)
    #     ax.set_ylim(0., 1.)
    #     line, = ax.plot(0., 0.)
    #
    #     x_data = []
    #     y_data = []
    #
    #     def animation_frame(i: int):
    #         """Отрисовка линии до i_end точки на j-ом слое"""
    #         x_data.append(self.time_list[j][i])
    #         y_data.append(self.flux_meas[j][i])
    #
    #         line.set_xdata(x_data)
    #         line.set_ydata(y_data)
    #         return line,
    #
    #     animation = FuncAnimation(fig, func=animation_frame, frames=range(x_len), interval=1)
    #     plt.show()