import pyomo.environ as pyomo
from pyomo.common.fileutils import this_file_dir
import sys, os, json
import pandas as pd
import numpy as np

# add rankine cycle code to path
sys.path.append(os.path.join(this_file_dir(),"../../models/simple_case"))
import simple_rankine_cycle as src


revenue_prescient=[]
nstartups_prescient=[]
capfactor_prescient=[]
operating_cost_prescient=[]

folder = 'pricetaker_2'
coal_price = 50.0
startup_cost = 61.09

with open('_verification_runs/{}/parameters.json'.format(folder), 'r') as infile:
    data = json.load(infile)


result_data_dir = os.path.join(this_file_dir(),"_verification_runs/{}/deterministic_simulation_output_index_0".format(folder))
bus_detail_df = pd.read_csv(os.path.join(result_data_dir,'bus_detail.csv'))
thermal_detail_df = pd.read_csv(os.path.join(result_data_dir,'thermal_detail.csv'))
coal_generator = '123_STEAM_3'
gen_results = thermal_detail_df.loc[thermal_detail_df['Generator'] == coal_generator] #output results for this generator
dispatch = gen_results["Dispatch"]
pmax=data['pmax']

# we don't have the results file, but this is what the capex is for price taker
cap_cost = 408.5693287301804

# This takes awhile...
op_costs = []
for output_value in dispatch.to_numpy():
    m = src.square_problem(
        heat_recovery=True,
        capital_fs=False,
        calc_boiler_eff=True,
        p_max=pmax,     # The pmax used for efficiency calculations
        net_power=output_value, # The operating output (used as nameplate for capex)
        coal_price=coal_price
    )
    op_costs.append(pyomo.value(m.fs.operating_cost))


# get LMPs
bus_results = bus_detail_df.loc[bus_detail_df['Bus'] == "CopperSheet"] #output results for this generator
lmp = np.copy(bus_results["LMP"])
lmp[lmp > 200] = 200
lmp_da = np.copy(bus_results["LMP DA"])

##############################################
# Calculate total revenue
##############################################
#more complex revenue, not really captures in Prescient runs
dispatch_da = gen_results["Dispatch DA"]
uplift_payment = gen_results['Unit Uplift Payment']
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
total_op_cost = sum(op_costs)/1e6
total_dispatch = dispatch.sum()
cap_factor = total_dispatch/(pmax*8736)

total_startup_cost = startup_cost*n_startups*pmax/1e6
number_of_years = 20
net_revenue = number_of_years*(total_revenue - total_op_cost - total_startup_cost) - cap_cost 

data_summary = {
    'total_revenue_yr':total_revenue,
    'total_op_cost_yr':total_op_cost,
    'total_dispatch_yr':total_dispatch,
    'cap_factor':cap_factor,
    'cap_cost':cap_cost,
    'startup_cost':total_startup_cost,
    'net_revenue':net_revenue
}

# write to JSON
with open('_verification_runs/{}_summary.json'.format(folder), 'w') as outfile:
    json.dump(data_summary, outfile)



# test out different revenue op cost combinations
total_revenue = 14.29 
total_op_cost = 9.0
total_startup_cost = 0.81
net_revenue = number_of_years*(total_revenue - total_op_cost - total_startup_cost) - cap_cost
print(net_revenue)


total_revenue = 12.3 
total_op_cost = 9.7
total_startup_cost = 0.71
net_revenue = number_of_years*(total_revenue - total_op_cost - total_startup_cost) - cap_cost
print(net_revenue)

total_revenue = 30.08
total_op_cost = 20.8
total_startup_cost = 0.01
net_revenue = number_of_years*(total_revenue - total_op_cost - total_startup_cost) - cap_cost
print(net_revenue)

total_revenue = 30.08
total_op_cost = 20.8
total_startup_cost = 0.01
net_revenue = number_of_years*(total_revenue - total_op_cost - total_startup_cost) - cap_cost
print(net_revenue)



total_revenue = 16.57
total_op_cost = 13.1
total_startup_cost = 0.0
net_revenue = number_of_years*(total_revenue - total_op_cost - total_startup_cost) - cap_cost
print(net_revenue)



# ALAMO surrogate
total_revenue = 26.6 
total_op_cost = 16.0
total_startup_cost = 0.81
net_revenue = number_of_years*(total_revenue - total_op_cost - total_startup_cost) - cap_cost
print(net_revenue)