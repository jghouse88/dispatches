import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib
import pandas as pd
matplotlib.rc('font', size=12)
plt.rc('axes', titlesize=12)
from pyomo.common.fileutils import this_file_dir
import matplotlib.ticker as ticker

data_dir = os.path.join(this_file_dir(),"../prescient_summary_dataframes")

df1 = pd.read_hdf(os.path.join(data_dir,"prescient_generator_inputs.h5"))
df2 = pd.read_hdf(os.path.join(data_dir,"prescient_generator_outputs.h5"))
df3 = pd.read_hdf(os.path.join(data_dir,"prescient_generator_startups.h5"))

p_max_vector = df1["p_max"].unique()

case1 = [0.3, 1.0, 4, 1, 15.0, 0.0, 0]
case2 = [0.3, 1.0, 4, 1, 15.0, 2.5, 0]
case3 = [0.3, 1.0, 4, 1, 15.0, 0.0, 4]
case4 = [0.3, 1.0, 4, 1, 15.0, 2.5, 4]
case5 = [0.3, 1.0, 4, 1, 25.0, 0.0, 0]
case6 = [0.3, 1.0, 4, 1, 25.0, 2.5, 0]
case7 = [0.3, 1.0, 4, 1, 25.0, 0.0, 4]
case8 = [0.3, 1.0, 4, 1, 25.0, 2.5, 4]

cases = [case1, case2, case3, case4, case5, case6, case7, case8]

parm_dicts = []
for gen_case in cases:
	pmax_revenue = []
	pmax_cap_factor = []
	pmax_nstartups = []
	for pmax in p_max_vector:
		# prepend pmax to case list
	    case = [pmax] + gen_case
	    df_case = df1.loc[(df1.p_max == case[0]) & (df1.p_min_multi == case[1]) & \
	    (df1.ramp_multi == case[2]) & (df1.min_up == case[3]) & (df1.min_dn_multi == case[4]) & \
	    (df1.marginal_cost == case[5]) & (df1.no_load_cost == case[6]) & \
	    (df1.startup_cost_profile == case[7])]

	    pmax_revenue.append((df2.iloc[df_case.index]["Total Revenue [$]"].to_numpy()/1e6)[0])
	    pmax_cap_factor.append((df2.iloc[df_case.index]["Total Dispatch [MWh]"].to_numpy()/(pmax*8736))[0])
	    pmax_nstartups.append((df3.iloc[df_case.index]["# Startups"].to_numpy())[0])

	parm_dict = {}
	parm_dict['cap_factor'] = pmax_cap_factor
	parm_dict['revenue'] = pmax_revenue
	parm_dict['n_startups'] = pmax_nstartups
	parm_dicts.append(parm_dict)


colors = matplotlib.cm.tab10(range(8))
markers = ['o','o','o','o','*','*','*','*']

#revenue
fig,axs = plt.subplots(1,3)
axs[0].set_xlabel("$P_{max}$ [MW]", fontweight='bold', fontsize=15)
axs[0].set_ylabel("Revenue [MM$]", fontweight='bold', fontsize=15)
axs[1].set_xlabel("$P_{max}$ [MW]", fontweight='bold', fontsize=15)
axs[1].set_ylabel("Capacity Factor", fontweight='bold', fontsize=15)
axs[2].set_xlabel("$P_{max}$ [MW]", fontweight='bold', fontsize=15)
axs[2].set_ylabel("# Startups", fontweight='bold', fontsize=15)

for p in range(len(parm_dicts)):
    revenue = parm_dicts[p]['revenue']
    axs[0].scatter(p_max_vector, revenue, color=colors[p], marker=markers[p], alpha=0.8, s=120)


    cap_fac = parm_dicts[p]['cap_factor']
    axs[1].scatter(p_max_vector, cap_fac, color=colors[p], marker=markers[p], alpha=0.8, s=120)

    n_startups = parm_dicts[p]['n_startups']
    axs[2].scatter(p_max_vector, n_startups, color=colors[p], marker=markers[p], alpha=0.8, s=120)


for ax in axs: 
	ratio = 0.8
	x_left, x_right = ax.get_xlim()
	y_low, y_high = ax.get_ylim()
	ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
	ax.tick_params(direction="in", top=True, right=True, labelsize=15)



plt.figlegend(['Case1', 'Case2', 'Case3', 'Case4', 'Case5', 'Case6', 'Case7', 'Case8'], loc=2, prop={'size': 15}, title="", ncol=4)
fig.set_size_inches(12,4) 
plt.tight_layout()
fig.subplots_adjust(left=0.05, bottom=0.15, right=0.97, top=0.80)  


plt.savefig("sweep_figures/scenario_sweeps.png")
plt.savefig("sweep_figures/scenario_sweeps.pdf")

