#Calculate correlations on unfiltered data
import pandas as pd
import os
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('font', size=24)
plt.rc('axes', titlesize=24)     # fontsize of the axes title
import numpy as np
import seaborn as sns
import scipy
# import pingouin as pg

f_perturbed_inputs = os.path.join(os.getcwd(),"../prescient_data/prescient_generator_inputs.h5")
f_perturbed_outputs = os.path.join(os.getcwd(),"../prescient_data/prescient_generator_outputs.h5")
f_startups = os.path.join(os.getcwd(),"../prescient_data/prescient_generator_startups.h5")

df_perturbed_inputs = pd.read_hdf(f_perturbed_inputs)
df_perturbed_outputs = pd.read_hdf(f_perturbed_outputs)
df_nstartups = pd.read_hdf(f_startups)

df_revenue = df_perturbed_outputs["Total Revenue [$]"] / 1e6
ax = df_revenue.plot.hist(bins=100, alpha=0.5)

plt.xlabel("Annual Revenue [MM$]")
plt.ylabel("")


plt.savefig("revenue_histogram.svg")

plt.show()
