# импорт классов
from .calc_flux import Wave
from .deposition_simulation import SetUpParameters
from .design_class import Design
from .target import Target
from .c_values import ProcessedStatistics
from .simulation_info import SimInfo
from .sim_params import SimParams
from .statistics_info import StatInfo

# импорт вычисляющих функций
from .calc_flux import calc_flux, grad_flux
from .deposition_simulation import simulation
from .statistics_info import mean_error_norm

# импорт функций строящих графики
from .correlation_plots import sigmas_plot
from .monitoring_curve_plot import monitoring_curve_plot
from .visualisation import thickness_error_box_plot, thickness_error_violin_plot, c_hist