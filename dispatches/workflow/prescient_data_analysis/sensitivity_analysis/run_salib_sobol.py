#Use highest accuracy surrogate with all terms
import pandas as pd
import os
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('font', size=18)
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('text', usetex=True)
import numpy as np
import pickle
import math
import json

from sklearn.model_selection import  train_test_split
from SALib.sample import saltelli         #saltelli sampler
from SALib.analyze import sobol
import SALib

f_perturbed_inputs = os.path.join(os.getcwd(),"../prescient_data/prescient_generator_inputs.h5")
f_perturbed_outputs = os.path.join(os.getcwd(),"../prescient_data/prescient_generator_outputs.h5")
f_startups = os.path.join(os.getcwd(),"../prescient_data/prescient_generator_startups.h5")

df_perturbed_inputs = pd.read_hdf(f_perturbed_inputs)
df_perturbed_outputs = pd.read_hdf(f_perturbed_outputs)
df_nstartups = pd.read_hdf(f_startups)

df_inputs = df_perturbed_inputs.iloc[:,[1,2,3,4,5,6,7,9]]
cols = df_inputs.columns.to_list()
df_bounds = [[min(df_inputs.iloc[:,i]),max(df_inputs.iloc[:,i])] for i in range(0,8)]

predicted_revenue = df_perturbed_outputs["Total Revenue [$]"]

########################
#REVENUE
########################
x = df_inputs.to_numpy()
z = predicted_revenue.to_numpy()/1e6
X_train, X_test, z_train, z_test = train_test_split(x, z, test_size=0.33, random_state=42)

xm = np.mean(X_train,axis = 0)
xstd = np.std(X_train,axis = 0)
zm = np.mean(z_train)
zstd = np.std(z_train)

x_scaled = (x - xm) / xstd
z_scaled = (z - zm) / zstd

with open("../surrogate_models/scikit/models/scikit_revenue.pkl", 'rb') as f:
    model = pickle.load(f)

#SALib problem dict
problem = {'num_vars': 8,'names': cols,'bounds': df_bounds}
param_values = saltelli.sample(problem, 2**14, calc_second_order=True)

scaled_param_values = (param_values - xm)/xstd
Y  = np.array([model.predict(inputs.reshape(1,8))[0] for inputs in scaled_param_values])
Y_unscaled = Y*zstd + zm

#Perform Sobol analysis
Si_sobol = sobol.analyze(problem, Y, calc_second_order=True)

S1 = Si_sobol["S1"]
S1_conf = Si_sobol["S1_conf"]
ST = Si_sobol["ST"]
ST_conf = Si_sobol["ST_conf"]
S2 = Si_sobol["S2"]
S2_conf = Si_sobol["S2_conf"]

# check sums
S2f = S2.flatten()
S2_values = [v for v in S2f if not np.isnan(v)]
print("Sum of Sobol indices (revenue): ", sum(S1) + sum(S2_values))

x = ["X{}".format(i) for i in range(1,9)]
Si_results = {"X":x,"S1":list(S1),"S1_conf":list(S1_conf),"ST":list(ST)}#,"S2":list(S2)}
filename = os.path.join("salib_results","Si_sobol_sensitivies_revenue.json")
with open(filename, 'w') as f:
    json.dump(Si_results, f)

S1T = np.array([S1, S1_conf, ST, ST_conf])
np.savetxt('salib_results/sobol_revenue_1.csv', S1T, delimiter='&',fmt='%.2e',newline="\\\ \n")

# write S2 results
x = []
y = []
e = []
for i in range(1,9):
    for j in range(i,9):
        if i != j:
            x.append(r"$X_{}X_{}$".format(i,j))
            y.append("{:.2e}".format(S2[i-1][j-1]))
            e.append("{:.2e}".format(Si_sobol["S2_conf"][i-1][j-1]))
S2_results = np.array(list(zip(x,y,e)))
np.savetxt('salib_results/sobol_revenue_2.csv', S2_results, delimiter='&', fmt='%s', newline="\\\ \n")


#Make plots (REVENUE)
#First order
x = ["X{}".format(i) for i in range(1,9)]
y = S1
e = Si_sobol["S1_conf"]
fig = plt.figure(figsize = (12,12))
fig.subplots_adjust(left = 0.2)
ax = fig.add_subplot(111)
ax.errorbar(x, y, yerr=e, fmt='o',markersize = 20,capsize = 10)
ax.set_ylabel("Sensitivity Index [Revenue]")
ax.set_xlabel("Input Parameter")
fig.savefig("salib_results/Si_sobol_1.png")

#Total order
y = ST
e = Si_sobol["ST_conf"]
fig = plt.figure(figsize = (12,12))
fig.subplots_adjust(left = 0.2)
ax = fig.add_subplot(111)
ax.errorbar(x, y, yerr=e, fmt='o',markersize = 20,capsize = 10)
ax.set_ylabel("Sensitivity Index [Revenue]")
ax.set_xlabel("Input Parameter")
fig.savefig("salib_results/Si_sobol_T.png")

#Second order
matplotlib.rc('font', size=12)
plt.rc('axes', titlesize=12)
x = []
y = []
e = []
for i in range(1,9):
    for j in range(i,9):
        if i != j:
            x.append("X{}X{}".format(i,j))
            y.append(S2[i-1][j-1])
            e.append(Si_sobol["S2_conf"][i-1][j-1])

fig = plt.figure(figsize = (24,8))
fig.subplots_adjust(left = 0.2)
ax = fig.add_subplot(111)
ax.errorbar(x, y, yerr=e, fmt='o',markersize = 20,capsize = 10)
ax.set_ylabel("Sensitivity Index \n [Revenue]")
ax.set_xlabel("Input Parameter")
plt.xticks(rotation = 45)
fig.tight_layout()
fig.savefig("salib_results//Si_sobol_2.png")


########################
#N-Startups
########################
x = df_inputs.to_numpy()
z = df_nstartups["# Startups"].to_numpy()
X_train, X_test, z_train, z_test = train_test_split(x, z, test_size=0.33, random_state=42)

xm = np.mean(X_train,axis = 0)
xstd = np.std(X_train,axis = 0)
zm = np.mean(z_train)
zstd = np.std(z_train)

x_scaled = (x - xm) / xstd
z_scaled = (z - zm) / zstd

with open("../surrogate_models/scikit/models/scikit_nstartups.pkl", 'rb') as f:
    model = pickle.load(f)

problem = {'num_vars': 8,'names': cols,'bounds': df_bounds}
param_values = saltelli.sample(problem, 2**14, calc_second_order=True)
scaled_param_values = (param_values - xm)/xstd
Y  = np.array([model.predict(inputs.reshape(1,8))[0] for inputs in scaled_param_values])
Y_unscaled = Y*zstd + zm

#Perform Sobol analysis
Si_sobol = sobol.analyze(problem, Y, calc_second_order=True)

S1 = Si_sobol["S1"]
S1_conf = Si_sobol["S1_conf"]
ST = Si_sobol["ST"]
ST_conf = Si_sobol["ST_conf"]
S2 = Si_sobol["S2"]
S2_conf = Si_sobol["S2_conf"]

# check sums
S2f = S2.flatten()
S2_values = [v for v in S2f if not np.isnan(v)]
print("Sum of Sobol indices (revenue): ", sum(S1) + sum(S2_values))

x = ["X{}".format(i) for i in range(1,9)]
Si_results = {"X":x,"S1":list(S1),"S1_conf":list(S1_conf),"ST":list(ST)}#,"S2":list(S2)}
filename = os.path.join("salib_results","Si_sobol_sensitivies_revenue.json")
with open(filename, 'w') as f:
    json.dump(Si_results, f)

S1T = np.array([S1, S1_conf, ST, ST_conf])
np.savetxt('salib_results/sobol_startups_1.csv', S1T, delimiter='&',fmt='%.2e',newline="\\\ \n")

# write S2 results
x2 = []
y2 = []
e2 = []
for i in range(1,9):
    for j in range(i,9):
        if i != j:
            x2.append(r"$X_{}X_{}$".format(i,j))
            y2.append("{:.2e}".format(S2[i-1][j-1]))
            e2.append("{:.2e}".format(Si_sobol["S2_conf"][i-1][j-1]))
S2_results = np.array(list(zip(x2,y2,e2)))
np.savetxt('salib_results/sobol_startups_2.csv', S2_results, delimiter='&', fmt='%s', newline="\\\ \n")

S2_results_all = np.array(list(zip(x2,y,e,y2,e2)))
np.savetxt('salib_results/sobol_rev_and_startups_2.csv', S2_results_all, delimiter='&', fmt='%s', newline="\\\ \n")

# #Perform Sobol analysis
# Si_sobol = sobol.analyze(problem, Y,calc_second_order = True)
#
# S1 = Si_sobol["S1"]
# S1_conf = Si_sobol["S1_conf"]
# ST = Si_sobol["ST"]
# S2 = Si_sobol["S2"]
#
# x = ["X{}".format(i) for i in range(1,9)]
# Si_results = {"X":x,"S1":list(S1),"S1_conf":list(S1_conf),"ST":list(ST)}#,"S2":list(S2)}
# filename = "salib_results/Si_sobol_sensitivies_nstartups.json"
# with open(filename, 'w') as f:
#     json.dump(Si_results, f)
#
# np.savetxt('salib_results/sobol_nstart.csv', ST.reshape(1, ST.shape[0]).round(2), delimiter='&',fmt='%.2g')
