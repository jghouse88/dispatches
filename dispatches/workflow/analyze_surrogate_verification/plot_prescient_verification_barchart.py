import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('font', size=12)
plt.rc('axes', titlesize=12)
from itertools import cycle, islice
from pyomo.common.fileutils import this_file_dir
import os
import json
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
##############################################
#Prescient Simulation Results
##############################################

dispatch_zones_surrogate=[]
dispatch_zones_prescient=[]

revenue_surrogate=[]
revenue_prescient=[]

nstartups_surrogate=[]
nstartups_prescient=[]

capfactor_surrogate=[]
capfactor_prescient=[]

zone_outputs = np.array([0.0,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,1.0])
folders = ['alamo_run_0','scikit_run_0','scikit_run_1','scikit_run_2']
for folder in folders:
    result_data_dir = os.path.join(this_file_dir(),"_verification_runs/{}/deterministic_simulation_output_index_0".format(folder))
    surrogate_data_filename = os.path.join(this_file_dir(),"_verification_runs/{}/verification_inputs.json".format(folder))

    with open(surrogate_data_filename) as f:
        surrogate_data = json.load(f)

    bus_detail_df = pd.read_csv(os.path.join(result_data_dir,'bus_detail.csv'))
    thermal_detail_df = pd.read_csv(os.path.join(result_data_dir,'thermal_detail.csv'))

    #Look at an example coal generator (This is the generator we studied for the big Prescient simulation data set)
    coal_generator = '123_STEAM_3'
    gen_results = thermal_detail_df.loc[thermal_detail_df['Generator'] == coal_generator] #output results for this generator
    dispatch = gen_results["Dispatch"]

    # bus_id = np.array(coal_gen_data["Bus ID"])[0]
    # bus_name = np.array(bus_df.loc[bus_df['Bus ID'] == bus_id]["Bus Name"])[0]
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
    #dam_revenue = np.nanprod(np.vstack((lmp_da,dispatch_da.values)),axis = 0)
    #revenue = np.nansum(np.vstack((rtm_revenue,dam_revenue,uplift_payment.values)),axis = 0)
    total_revenue = sum(rtm_revenue)/1e6
    total_dispatch = dispatch.sum()


    ##############################################
    # Calculate dispatch zones
    ##############################################
    pmax = surrogate_data["market_inputs"][0]
    pmin = surrogate_data["market_inputs"][1]*pmax

    bins_col_names = ["# Hours 0%","# Hours 0-10%","# Hours 10-20%","# Hours 20-30%","# Hours 30-40%","# Hours 40-50%","# Hours 50-60%","# Hours 60-70%","# Hours 70-80%","# Hours 80-90%","# Hours 90-100%"]
    dispatch_ranges = [(0,10),(10,20),(20,30),(30,40),(40,50),(50,60),(60,70),(70,80),(80,90),(90,100)]

    prescient_zone_hours = np.zeros(len(dispatch_ranges)+1)
    dispatch_prescient = dispatch.to_numpy()
    n_hours_off = len(np.nonzero(dispatch_prescient == 0)[0])

    dispatch_on = dispatch_prescient[dispatch_prescient > 0.0]
    scaled_dispatch = np.round((dispatch_on - pmin) / (pmax - pmin) * 100,3)

    assert(n_hours_off + len(scaled_dispatch) == 8736)
    prescient_zone_hours[0] = n_hours_off
    for i in range(len(dispatch_ranges)-1):
        zone = dispatch_ranges[i]
        n_hours_zone = len(np.where((scaled_dispatch >= zone[0]) & (scaled_dispatch < zone[1]))[0])
        prescient_zone_hours[i+1] = n_hours_zone
    prescient_zone_hours[-1] = len(np.where((scaled_dispatch >= dispatch_ranges[-1][0]) & (scaled_dispatch <= dispatch_ranges[-1][1]))[0])

    assert(sum(prescient_zone_hours) == 8736)



    ##############################################
    # Calculate number of startups
    ##############################################
    thermal_detail_diff = gen_results['Unit State'].diff().fillna(0)
    state_changes = thermal_detail_diff.value_counts()
    if 1 in state_changes.keys():
        n_startups = state_changes[1]
    else:
        n_startups = 0

    dispatch_zones_surrogate.append(surrogate_data['dispatch_zones'])
    dispatch_zones_prescient.append(prescient_zone_hours)

    revenue_surrogate.append(surrogate_data['revenue_surrogate'])
    revenue_prescient.append(total_revenue)

    nstartups_surrogate.append(surrogate_data['nstartups_per_year'])
    nstartups_prescient.append(n_startups)

    raw_zone_outputs = zone_outputs*(pmax-pmin) + pmin
    s_dispatch_zones = surrogate_data['dispatch_zones']
    zone_hours_on = s_dispatch_zones[1:]
    capfactor_surrogate.append((np.dot(zone_hours_on,raw_zone_outputs) / (pmax*8736)))
    capfactor_prescient.append(total_dispatch/(pmax*8736))

# Use Dataframe to get the plot
# fig, ax = plt.subplots(figsize = (16,8))
fig, ax = plt.subplots()
ax.set_ylabel("Hours at Power Output", fontweight='bold', fontsize=18)
ax.set_xlabel("Scaled Power Output (% of maximum)", fontweight='bold', fontsize=18)
labels = ["Off","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"]

fig.set_size_inches(12,9)

columns = list(zip([labels[0],labels[1],labels[10]],
    list([dispatch_zones_surrogate[0][0],dispatch_zones_surrogate[0][1],dispatch_zones_surrogate[0][10]]),
    list([dispatch_zones_prescient[0][0],dispatch_zones_prescient[0][1],dispatch_zones_prescient[0][10]]),
    list([dispatch_zones_surrogate[1][0],dispatch_zones_surrogate[1][1],dispatch_zones_surrogate[1][10]]),
    list([dispatch_zones_prescient[1][0],dispatch_zones_prescient[1][1],dispatch_zones_prescient[1][10]]),
    list([dispatch_zones_surrogate[2][0],dispatch_zones_surrogate[2][1],dispatch_zones_surrogate[2][10]]),
    list([dispatch_zones_prescient[2][0],dispatch_zones_prescient[2][1],dispatch_zones_prescient[2][10]]),
    list([dispatch_zones_surrogate[3][0],dispatch_zones_surrogate[3][1],dispatch_zones_surrogate[3][10]]),
    list([dispatch_zones_prescient[3][0],dispatch_zones_prescient[3][1],dispatch_zones_prescient[3][10]]),
    ))

df_plt = pd.DataFrame(columns, columns=["Scaled Power Output (% of maximum)",
    "ALAMO-1","Verify ALAMO-1","NN-1","Verify NN-1","NN-2","Verify NN-2","NN-3","Verify NN-3"])
my_colors = list(islice(cycle(['tab:blue','tab:blue','tab:green' ,'tab:green','tab:purple','tab:purple','tab:red','tab:red']), None, len(df_plt)*8))
df_plt.plot(ax=ax, x="Scaled Power Output (% of maximum)", kind='bar', stacked=False, color=my_colors, legend=False, width = 0.9)
ax.tick_params(axis='x', labelrotation = 0)


#hacky way to set bar hatches
bars = ax.patches
pattern = "//"
n_groups = 3
hatch_indices = np.array([np.arange(n_groups,n_groups+n_groups,1)+n_groups*i*2 for i in range(4)]).flatten()
for bar in np.array(bars)[hatch_indices]:
    bar.set_hatch(pattern)

legend_elements = [Patch(facecolor='tab:blue', edgecolor='tab:blue', label='ALAMO-1')]
legend_elements.append(Patch(facecolor='tab:green', edgecolor='tab:green', label='NN-1'))
legend_elements.append(Patch(facecolor='tab:purple', edgecolor='tab:purple', label='NN-3'))
legend_elements.append(Patch(facecolor='tab:red', edgecolor='tab:red', label='NN-2'))
legend_elements.append(Line2D([],[],linestyle=''))
legend_elements.append(Patch(facecolor='white', edgecolor='black', label='Surrogate Solution'))
legend_elements.append(Patch(facecolor='white', edgecolor='black', hatch = "//", label='Prescient Verification'))
ax.legend(handles=legend_elements, loc='upper left', prop={'size': 18})
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.tight_layout()
plt.savefig("zone_verification_all.png")
plt.savefig("zone_verification_all.pdf")



 # plot intermediate zones
fig, ax = plt.subplots()
ax.set_ylabel("Hours at Power Output", fontweight='bold', fontsize=18)
ax.set_xlabel("Scaled Power Output (% of maximum)", fontweight='bold', fontsize=18)
labels = ["Off","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"]
fig.set_size_inches(12,9)

columns = list(zip(labels[2:10],
    list(dispatch_zones_surrogate[0][2:10]),list(dispatch_zones_prescient[0][2:10]),
    list(dispatch_zones_surrogate[1][2:10]),list(dispatch_zones_prescient[1][2:10]),
    list(dispatch_zones_surrogate[2][2:10]),list(dispatch_zones_prescient[2][2:10]),
    list(dispatch_zones_surrogate[3][2:10]),list(dispatch_zones_prescient[3][2:10])
    ))

df_plt = pd.DataFrame(columns, columns=["Scaled Power Output (% of maximum)",
    "ALAMO-1","Verify ALAMO-1","NN-1","Verify NN-1","NN-2","Verify NN-2","NN-3","Verify NN-3"])
my_colors = list(islice(cycle(['tab:blue','tab:blue','tab:green' ,'tab:green','tab:purple','tab:purple','tab:red','tab:red']), None, len(df_plt)*8))
df_plt.plot(ax = ax, x="Scaled Power Output (% of maximum)", kind='bar', stacked=False, color=my_colors, legend=False, width = 0.9)
ax.tick_params(axis='x', labelrotation = 0)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

#hacky way to set bar hatches
bars = ax.patches
pattern = "//"
n_groups = 8
hatch_indices = np.array([np.arange(n_groups,n_groups+n_groups,1)+n_groups*i*2 for i in range(4)]).flatten()
for bar in np.array(bars)[hatch_indices]:
    bar.set_hatch(pattern)

legend_elements = [Patch(facecolor='tab:blue', edgecolor='tab:blue', label='ALAMO-1')]
legend_elements.append(Patch(facecolor='tab:green', edgecolor='tab:green', label='NN-1'))
legend_elements.append(Patch(facecolor='tab:purple', edgecolor='tab:purple', label='NN-3'))
legend_elements.append(Patch(facecolor='tab:red', edgecolor='tab:red', label='NN-2'))
legend_elements.append(Line2D([],[],linestyle=''))
legend_elements.append(Patch(facecolor='white', edgecolor='black', label='Surrogate Solution'))
legend_elements.append(Patch(facecolor='white', edgecolor='black', hatch = "//", label='Prescient Verification'))
ax.legend(handles=legend_elements, loc='upper left', prop={'size': 18})
plt.tight_layout()

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)


plt.savefig("zone_verification_mid.png")
plt.savefig("zone_verification_mid.pdf")
