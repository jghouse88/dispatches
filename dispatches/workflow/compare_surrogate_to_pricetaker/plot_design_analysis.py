# this script creates the key design figures that show why surrogates go beyond price taker

import sys, os
sys.path.append("../")
import pandas as pd

import numpy as np
from copy import copy
from get_surrogate_results import predict_with_surrogate
from pyomo.common.fileutils import this_file_dir

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib
matplotlib.rc('font', size=12)
plt.rc('axes', titlesize=12)

data_dir = os.path.join(this_file_dir(),"../prescient_summary_dataframes")
df1 = pd.read_hdf(os.path.join(data_dir,"prescient_generator_inputs.h5"))
df2 = pd.read_hdf(os.path.join(data_dir,"prescient_generator_outputs.h5"))

startup_csts = [0., 49.66991167, 61.09068702, 101.4374234,  135.2230393]
surrogate_mrg_csts = np.linspace(5,30,100)

#load price taker (vanilla) results pmax=177.5
df_price_taker_free = pd.read_pickle('_rminlp_price_taker_revenue_vs_marginal_cost_optimal_p_max.pkl')
mrg_csts_price_taker_free = df_price_taker_free["marginal_cost ($/MWh)"].to_numpy()
rev_price_taker_free = df_price_taker_free["annual_revenue ($M/yr)"].to_numpy()
cap_price_taker_free = df_price_taker_free["avg_capacity_factor (%)"].to_numpy()/100

#load price taker (vanilla) results pmax=177.5
df_price_taker_low = pd.read_pickle('_rminlp_price_taker_revenue_vs_marginal_cost_fixed_p_max_177.pkl')
mrg_csts_price_taker_low = df_price_taker_low["marginal_cost ($/MWh)"].to_numpy()
rev_price_taker_low = df_price_taker_low["annual_revenue ($M/yr)"].to_numpy()
cap_price_taker_low = df_price_taker_low["avg_capacity_factor (%)"].to_numpy()/100

#pmax=266.75
df_price_taker_mid = pd.read_pickle('_rminlp_price_taker_revenue_vs_marginal_cost_fixed_p_max_266.pkl')
mrg_csts_price_taker_mid = df_price_taker_mid["marginal_cost ($/MWh)"].to_numpy()
rev_price_taker_mid = df_price_taker_mid["annual_revenue ($M/yr)"].to_numpy()
cap_price_taker_mid = df_price_taker_mid["avg_capacity_factor (%)"].to_numpy()/100

#pmax=443.75
df_price_taker_high = pd.read_pickle('_rminlp_price_taker_revenue_vs_marginal_cost_pmax_443.pkl')
mrg_csts_price_taker_high = df_price_taker_high["marginal_cost ($/MWh)"].to_numpy()
rev_price_taker_high = df_price_taker_high["annual_revenue ($M/yr)"].to_numpy()
cap_price_taker_high = df_price_taker_high["avg_capacity_factor (%)"].to_numpy()/100

#CASE 1: low pmax
cases1 = []
surrogate_data1 = []
pmax = 177.5
for i in range(5):
    start_cst_index=i
    case = [pmax, 0.3, 1.0, 4, 1, 0.0, start_cst_index, surrogate_mrg_csts]

    df_case = copy(df1.loc[(df1.p_max == case[0]) & (df1.p_min_multi == case[1]) & (df1.ramp_multi == case[2]) & \
    (df1.min_up == case[3]) & (df1.min_dn_multi == case[4]) & (df1.no_load_cost == case[5]) & (df1.startup_cost_profile == case[6])])

    df_case["Avg LMP"] = df2["Average LMP [$/MWh]"].to_numpy()[df_case.index]

    rev = df2.iloc[df_case.index]["Total Revenue [$]"].to_numpy()/1e6
    df_case["Revenue [$MM]"] = rev

    cap_factor = df2.iloc[df_case.index]["Total Dispatch [MWh]"].to_numpy()/(pmax*8736)
    df_case["Capacity Factor"] = cap_factor

    cases1.append(df_case)
    surrogate_case = copy(case)
    surrogate_case[6] = startup_csts[i]
    surrogate_data1.append(predict_with_surrogate(*surrogate_case))

sim_indices1 = [cases1[i].index.to_list() for i in range(len(cases1))]

#CASE 2: high pmax
cases2 = []
surrogate_data2 = []
pmax = 443.75
for i in range(5):
    start_cst_index=i
    case = [pmax, 0.3, 1.0, 4, 1, 0.0, start_cst_index, surrogate_mrg_csts]

    df_case = copy(df1.loc[(df1.p_max == case[0]) & (df1.p_min_multi == case[1]) & (df1.ramp_multi == case[2]) & \
    (df1.min_up == case[3]) & (df1.min_dn_multi == case[4]) & (df1.no_load_cost == case[5]) & (df1.startup_cost_profile == case[6])])
    
    df_case["Avg LMP"] = df2["Average LMP [$/MWh]"].to_numpy()[df_case.index]

    rev = df2.iloc[df_case.index]["Total Revenue [$]"].to_numpy()/1e6
    df_case["Revenue [$MM]"] = rev

    cap_factor = df2.iloc[df_case.index]["Total Dispatch [MWh]"].to_numpy()/(pmax*8736)
    df_case["Capacity Factor"] = cap_factor

    cases2.append(df_case)
    surrogate_case = copy(case)
    surrogate_case[6] = startup_csts[i]
    surrogate_data2.append(predict_with_surrogate(*surrogate_case))
sim_indices2 = [cases2[i].index.to_list() for i in range(len(cases2))]

#CASE 3: mid pmax
cases3 = []
surrogate_data3 = []
pmax = 266.25
for i in range(5):
    start_cst_index=i
    case = [pmax,0.3,1.0,4,1,0.0,start_cst_index,surrogate_mrg_csts]
    
    df_case = copy(df1.loc[(df1.p_max == case[0]) & (df1.p_min_multi == case[1]) & (df1.ramp_multi == case[2]) & \
    (df1.min_up == case[3]) & (df1.min_dn_multi == case[4]) & (df1.no_load_cost == case[5]) & (df1.startup_cost_profile == case[6])])
    
    df_case["Avg LMP"] = df2["Average LMP [$/MWh]"].to_numpy()[df_case.index]

    rev = df2.iloc[df_case.index]["Total Revenue [$]"].to_numpy()/1e6
    df_case["Revenue [$MM]"] = rev

    cap_factor = df2.iloc[df_case.index]["Total Dispatch [MWh]"].to_numpy()/(pmax*8736)
    df_case["Capacity Factor"] = cap_factor

    cases3.append(df_case)

    surrogate_case = copy(case)
    surrogate_case[6] = startup_csts[i]
    surrogate_data3.append(predict_with_surrogate(*surrogate_case))
sim_indices3 = [cases1[i].index.to_list() for i in range(len(cases1))]

revenues_pmax_low = [df["Revenue [$MM]"].to_numpy() for df in cases1]
revenues_pmax_surrogate_low = [surrogate_data1[i]["revenue"] for i in range(5)]
cap_factors_pmax_low = [df["Capacity Factor"].to_numpy() for df in cases1]
cap_factors_pmax_surrogate_low = [surrogate_data1[i]["capacity_factor"] for i in range(5)]

revenues_pmax_high = [df["Revenue [$MM]"].to_numpy() for df in cases2]
revenues_pmax_surrogate_high = [surrogate_data2[i]["revenue"] for i in range(5)]
cap_factors_pmax_high = [df["Capacity Factor"].to_numpy() for df in cases2]
cap_factors_pmax_surrogate_high = [surrogate_data2[i]["capacity_factor"] for i in range(5)]

revenues_pmax_mid = [df["Revenue [$MM]"].to_numpy() for df in cases3]
revenues_pmax_surrogate_mid = [surrogate_data3[i]["revenue"] for i in range(5)]
cap_factors_pmax_mid = [df["Capacity Factor"].to_numpy() for df in cases3]
cap_factors_pmax_surrogate_mid = [surrogate_data3[i]["capacity_factor"] for i in range(5)]

mrg_csts = [5,10,15,20,25,30]
colors = ["tab:purple","tab:blue","tab:green","tab:orange","tab:red"]


fig,axs = plt.subplots(figsize = (6,6))
axs.set_xlabel("Marginal Cost [$/MWh]", fontweight='bold', fontsize=18)
axs.set_ylabel("Revenue [$MM]", fontweight='bold', fontsize=18)
for i in range(5):
    plt.scatter(mrg_csts,revenues_pmax_low[i], s = 250, alpha = 0.75, color = colors[i])
    plt.plot(surrogate_mrg_csts,revenues_pmax_surrogate_low[i], alpha = 0.75, color = colors[i], linewidth = 2.0)
plt.scatter(mrg_csts_price_taker_low,rev_price_taker_low, s = 500, alpha = 1.0, marker = "X", color="black")
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[i], markersize=24, label='Start Profile {}'.format(i)) for i in range(5)]
legend_elements.append(Line2D([],[],linestyle=''))
legend_elements.append(Line2D([0], [0], marker='o', color="w", markerfacecolor="black" , markersize=24, label="Prescient Simulation"))
legend_elements.append(Line2D([0], [0], color = "black", linewidth=2.0, label = "NN Surrogate"))
legend_elements.append(Line2D([0], [0], marker='X', color="w", markerfacecolor="black", markersize=24,label="Price Taker"))
axs.legend(handles=legend_elements, loc='lower left', prop={'size': 15})
axs.tick_params(direction="in", top=True, right=True, labelsize=15)
plt.tight_layout()
plt.savefig("figures/revenue_mrg_cst_sweep_pmax_low.png")
plt.savefig("figures/revenue_mrg_cst_sweep_pmax_low.pdf")


fig,axs = plt.subplots(figsize = (6,6))
axs.set_xlabel("Marginal Cost [$/MWh]", fontweight='bold', fontsize=18)
axs.set_ylabel("Capacity Factor", fontweight='bold', fontsize=18)
for i in range(5):
    plt.scatter(mrg_csts,cap_factors_pmax_low[i], s = 250, alpha = 0.75, color = colors[i])
    plt.plot(surrogate_mrg_csts,cap_factors_pmax_surrogate_low[i], alpha = 0.75, color = colors[i], linewidth=2.0)
plt.scatter(mrg_csts_price_taker_low,cap_price_taker_low,label = "Price Taker", s = 500, alpha = 1.0, marker = "X", color="black")
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[i], markersize=24, label='Start Profile {}'.format(i)) for i in range(5)]
legend_elements.append(Line2D([],[],linestyle=''))
legend_elements.append(Line2D([0], [0], marker='o', color="w", markerfacecolor="black" , markersize=24, label="Prescient Simulation"))
legend_elements.append(Line2D([0], [0], color = "black", linewidth=2.0, label = "NN Surrogate"))
legend_elements.append(Line2D([0], [0], marker='X', color="w", markerfacecolor="black", markersize=24,label="Price Taker"))
axs.legend(handles=legend_elements, loc='lower left', prop={'size': 15})
axs.tick_params(direction="in", top=True, right=True, labelsize=15)
plt.tight_layout()
plt.savefig("figures/capfactor_mrg_cst_sweep_pmax_low.png")
plt.savefig("figures/capfactor_mrg_cst_sweep_pmax_low.pdf")


fig,axs = plt.subplots(figsize = (6,6))
axs.set_xlabel("Marginal Cost [$/MWh]", fontweight='bold', fontsize=18)
axs.set_ylabel("Revenue [$MM]", fontweight='bold', fontsize=18)
for i in range(5):
    plt.scatter(mrg_csts,revenues_pmax_high[i], s = 250, alpha = 0.75, color = colors[i])
    plt.plot(surrogate_mrg_csts,revenues_pmax_surrogate_high[i], alpha = 0.75, color = colors[i], linewidth=2.0)
plt.scatter(mrg_csts_price_taker_high,rev_price_taker_high,label = "Price Taker", s = 500, alpha = 1.0, marker = "X", color="black")
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[i], markersize=24, label='Start Profile {}'.format(i)) for i in range(5)]
legend_elements.append(Line2D([],[],linestyle=''))
legend_elements.append(Line2D([0], [0], marker='o', color="w", markerfacecolor="black" , markersize=24, label="Prescient Simulation"))
legend_elements.append(Line2D([0], [0], color = "black", linewidth=2.0, label = "NN Surrogate"))
legend_elements.append(Line2D([0], [0], marker='X', color="w", markerfacecolor="black", markersize=24,label="Price Taker"))
axs.legend(handles=legend_elements, loc='lower left', prop={'size': 15})
axs.tick_params(direction="in", top=True, right=True, labelsize=15)
plt.tight_layout()
plt.savefig("figures/revenue_mrg_cst_sweep_pmax_high.png")
plt.savefig("figures/revenue_mrg_cst_sweep_pmax_high.pdf")


fig,axs = plt.subplots(figsize = (6,6))
axs.set_xlabel("Marginal Cost [$/MWh]", fontweight='bold', fontsize=18)
axs.set_ylabel("Capacity Factor", fontweight='bold', fontsize=18)
for i in range(5):
    plt.scatter(mrg_csts,cap_factors_pmax_high[i], s = 250, alpha = 0.75, color = colors[i])
    plt.plot(surrogate_mrg_csts,cap_factors_pmax_surrogate_high[i], alpha = 0.75, color = colors[i], linewidth=2.0)
plt.scatter(mrg_csts_price_taker_high,cap_price_taker_high,label = "Price Taker", s = 500, alpha = 1.0, marker = "X", color="black")
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[i], markersize=24, label='Start Profile {}'.format(i)) for i in range(5)]
legend_elements.append(Line2D([],[],linestyle=''))
legend_elements.append(Line2D([0], [0], marker='o', color="w", markerfacecolor="black" , markersize=24, label="Prescient Simulation"))
legend_elements.append(Line2D([0], [0], color = "black", linewidth=2.0, label = "NN Surrogate"))
legend_elements.append(Line2D([0], [0], marker='X', color="w", markerfacecolor="black", markersize=24,label="Price Taker"))
axs.legend(handles=legend_elements, loc='lower left', prop={'size': 15})
axs.tick_params(direction="in", top=True, right=True, labelsize=15)
plt.tight_layout()
plt.savefig("figures/capfactor_mrg_cst_sweep_pmax_high.png")
plt.savefig("figures/capfactor_mrg_cst_sweep_pmax_high.pdf")


fig,axs = plt.subplots(figsize = (6,6))
axs.set_xlabel("Marginal Cost [$/MWh]", fontweight='bold', fontsize=18)
axs.set_ylabel("Revenue [$MM]", fontweight='bold', fontsize=18)
for i in range(5):
    plt.scatter(mrg_csts,revenues_pmax_mid[i], s = 250, alpha = 0.75, color = colors[i])
    plt.plot(surrogate_mrg_csts,revenues_pmax_surrogate_mid[i], alpha = 0.75, color = colors[i], linewidth=2.0)
plt.scatter(mrg_csts_price_taker_mid,rev_price_taker_mid,label = "Price Taker", s = 500, alpha = 1.0, marker = "X", color="black")
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[i], markersize=24, label='Start Profile {}'.format(i)) for i in range(5)]
legend_elements.append(Line2D([],[],linestyle=''))
legend_elements.append(Line2D([0], [0], marker='o', color="w", markerfacecolor="black" , markersize=24, label="Prescient Simulation"))
legend_elements.append(Line2D([0], [0], color = "black", linewidth=2.0, label = "NN Surrogate"))
legend_elements.append(Line2D([0], [0], marker='X', color="w", markerfacecolor="black", markersize=24,label="Price Taker"))
axs.legend(handles=legend_elements, loc='lower left', prop={'size': 15})
axs.tick_params(direction="in", top=True, right=True, labelsize=15)
plt.tight_layout()
plt.savefig("figures/revenue_mrg_cst_sweep_pmax_mid.png")
plt.savefig("figures/revenue_mrg_cst_sweep_pmax_mid.pdf")


fig,axs = plt.subplots(figsize = (6,6))
axs.set_xlabel("Marginal Cost [$/MWh]", fontweight='bold', fontsize=18)
axs.set_ylabel("Capacity Factor", fontweight='bold', fontsize=18)
for i in range(5):
    plt.scatter(mrg_csts,cap_factors_pmax_mid[i], s = 250, alpha = 0.75, color = colors[i])
    plt.plot(surrogate_mrg_csts,cap_factors_pmax_surrogate_mid[i], alpha = 0.75, color = colors[i], linewidth=2.0)
plt.scatter(mrg_csts_price_taker_mid,cap_price_taker_mid,label = "Price Taker", s = 500, alpha = 1.0, marker = "X", color="black")
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[i], markersize=24, label='Start Profile {}'.format(i)) for i in range(5)]
legend_elements.append(Line2D([],[],linestyle=''))
legend_elements.append(Line2D([0], [0], marker='o', color="w", markerfacecolor="black" , markersize=24, label="Prescient Simulation"))
legend_elements.append(Line2D([0], [0], color = "black", linewidth=2.0, label = "NN Surrogate"))
legend_elements.append(Line2D([0], [0], marker='X', color="w", markerfacecolor="black", markersize=24,label="Price Taker"))
axs.legend(handles=legend_elements, loc='lower left', prop={'size': 15})
axs.tick_params(direction="in", top=True, right=True, labelsize=15)
plt.tight_layout()
plt.savefig("figures/capfactor_mrg_cst_sweep_pmax_mid.png")
plt.savefig("figures/capfactor_mrg_cst_sweep_pmax_mid.pdf")
