# produce plot
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('font', size=18)
plt.rc('axes', titlesize=18)
import matplotlib.ticker as ticker
import os
import pickle
import json
import pandas as pd
from sklearn.model_selection import  train_test_split


f_inputs = os.path.join(os.getcwd(),"../../prescient_data/prescient_generator_inputs.h5")
f_outputs = os.path.join(os.getcwd(),"../../prescient_data/prescient_generator_outputs.h5")
df_inputs = pd.read_hdf(f_inputs)
df_outputs = pd.read_hdf(f_outputs)
predicted_revenue = df_outputs["Total Revenue [$]"]

x = df_inputs.iloc[:,[1,2,3,4,5,6,7,9]].to_numpy()
z = predicted_revenue.to_numpy()/1e6
X_train, X_test, z_train, z_test = train_test_split(x, z, test_size=0.33, random_state=42)

with open('models/training_parameters_revenue.json', 'r') as outfile:
    data = json.load(outfile)

xm = data['xm_inputs']
xstd = data['xstd_inputs']
zm = data['zm_revenue']
zstd = data['zstd_revenue']

#load up revenue model
with open("models/scikit_revenue.pkl", 'rb') as f:
    model = pickle.load(f)

with open('models/scikit_revenue_accuracy.json', 'r') as outfile:
    accuracy_dict = json.load(outfile)

X_test_scaled = (X_test - xm) / xstd
R2 = round(accuracy_dict["R2"],3)
predicted_revenue = model.predict(X_test_scaled)
predict_unscaled = predicted_revenue*zstd + zm

# plot results
plt.figure(figsize=(5,5))
plt.scatter(z_test, predict_unscaled, color="tab:blue", alpha = 0.01)
plt.plot([min(z), max(z)],[min(z), max(z)], color="black", linewidth=3.0)
plt.xlabel("True Revenue [MM$]", fontweight='bold')
plt.ylabel("Predicted Revenue [MM$]", fontweight='bold')
y_text = 0.75*(max(z) + min(z)) - min(z)
plt.annotate("$R^2 = {}$".format(R2),(0,y_text))


plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax = plt.gca()
ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
plt.tick_params(direction="in",top=True, right=True)

plt.tight_layout()
plt.savefig("figures/revenue_scikit.png")
plt.savefig("figures/revenue_scikit.pdf")
