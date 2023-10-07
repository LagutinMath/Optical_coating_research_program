# импорт классов
from .units import Wave, WidthForm, TermCase
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
from .process_statistics import process_statistics
from .tools import timer

# импорт функций строящих графики
from .monitoring_curve_plot import MonochromStrategyData
from .visualisation import thickness_error_box_plot, thickness_error_violin_plot, c_hist, sigmas_plot, std_values