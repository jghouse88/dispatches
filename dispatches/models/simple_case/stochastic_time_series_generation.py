# Import Pyomo libraries
from pyomo.environ import value
# from pyomo.util.infeasible import log_close_to_bounds

# from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util import get_solver
from pyparsing import lineStart

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from time import perf_counter


# RTS-GMLC Signal Unfiltered (new data - 10/10/2021)

with open('rts_results_all_prices_base_case.npy', 'rb') as f:
    dispatch = np.load(f)
    price = np.load(f)


day_1_clean = price[0:23]

# define number of scenarios for the noise in data
n_scenarios = 5

# noise signal
hours = {}

for i in range(0, len(day_1_clean)+1):
    hours[i] = np.random.normal(i, 5, n_scenarios)

df = pd.DataFrame(data=hours)
# plot the scenarios
# plt.plot(day_1_clean, c="r")
# for i in day_1_clean:
#     plt.plot(scenarios[i], linestyle="dashed")

# plt.show()

print(df)
