import pandas as pd
import numpy as np
from pyomo.common.fileutils import this_file_dir
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('font', size=12)
plt.rc('axes', titlesize=12)
import os


basecase_run = os.path.join(this_file_dir(),"basecase/_basecase_runs/_run_1/deterministic_simulation_output_index_0")

bus_df = pd.read_csv(os.path.join(basecase_run,"bus_detail.csv"))
gen_df = pd.read_csv(os.path.join(basecase_run,"thermal_detail.csv"))

coal_generator = '123_STEAM_3'
gen_results = gen_df.loc[gen_df['Generator'] == coal_generator]
dispatch = gen_results["Dispatch"]

pmax_dispatch = gen_results[dispatch == 350]
pmin_dispatch = gen_results[dispatch == 140]

op_cost_pmax = pmax_dispatch["Unit Cost"] / pmax_dispatch["Dispatch"]
op_cost_pmin = pmin_dispatch["Unit Cost"] / pmin_dispatch["Dispatch"]


# the average marginal cost from the basecase generator
mrg_cst = (op_cost_pmax.mode() + op_cost_pmin.mode()) / 2


bus_results = bus_df.loc[bus_df['Bus'] == "CopperSheet"] #output results for this generator
lmp = np.copy(bus_results["LMP"])
lmp[lmp > 200] = 200

pmax = 350
##############################################
# Calculate total revenue
##############################################
dispatch_da = gen_results["Dispatch DA"]
dispatch_diff = np.nansum(np.vstack((dispatch.values,-dispatch_da.values)),axis = 0)
rtm_revenue = np.nanprod(np.vstack((lmp,dispatch_diff)),axis = 0)

##############################################
# Calculate number of startups
##############################################
thermal_detail_diff = gen_results['Unit State'].diff().fillna(0)
state_changes = thermal_detail_diff.value_counts()
if 1 in state_changes.keys():
    n_startups = state_changes[1]
else:
    n_startups = 0

# calculate total values
total_revenue = sum(rtm_revenue)/1e6
total_dispatch = dispatch.sum()
cap_factor = total_dispatch/(pmax*8736)