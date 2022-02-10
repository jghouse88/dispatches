"""
Script to generate and save figures from results saved as a pandas dataframe
"""

__author__ = "Jaffer Ghouse"

from matplotlib import markers, pyplot as plt
import numpy as np
import pandas as pd

# load results saved as pandas dataframe
results = pd.read_pickle("price_taker_full_year_with_boil_eff_coal_50.pkl")
price=results["price"].tolist()
hour_list=results["hour_list"].tolist()
p_scenario=results["p_scenario"].tolist()
capacity_factor_dispatch=results["capacity_factor_dispatch"].tolist()

# Generate sub-plots: 4 sub-plots
# subplot(1,1): LMP and op cost vs. time
# subplot(1,2): LMP histogram
# subplot(2,1): Power schedule vs. time
# subplot(2,2): power profile histogram
"""
fig, axs = plt.subplots(2, 2)
# subplot(0,1): LMP vs. time
# axs[0].step(hour_list, price, linestyle="dotted", alpha=0.5, linewidth=0.75,
#                color="green")
axs[0, 0].scatter(hour_list, price, s=6, marker="o", facecolor="green",
                  linewidths=0.5, alpha=0.4, color="green")
axs[0, 0].grid(which="major", axis="both", linestyle="--")
axs[0, 0].set_xlabel("Time (h)")
axs[0, 0].set_ylabel("LMP ($/MWh)")

# subplot(1,0): p_scenario vs. time
# axs[1].step(hour_list, p_scenario, linestyle="dotted", alpha=0.5,
#             linewidth=0.75, color="blue")
axs[1, 0].scatter(hour_list, p_scenario, s=6, marker="o",
                  facecolor="blue", linewidths=0.5, alpha=0.4, color="blue")
axs[1, 0].grid(which="major", axis="both", linestyle="--")
axs[1, 0].set_xlabel("Time (h)")
axs[1, 0].set_ylabel("Power (MW)")

# subplot(0,1): price histogram
axs[0, 1].hist(price, bins=50, color="green")
axs[0, 1].grid(which="major", axis="both", linestyle="--")
axs[0, 1].set_xlabel("LMP ($/MWh)")
axs[0, 1].set_ylabel("Frequency")

# subplot(1,1): power histogram
# axs[1, 1].hist(p_scenario, color="blue")
# axs[1, 1].grid(which="major", axis="both", linestyle="--")
# axs[1, 1].set_xlabel("Power (MW)")
# axs[1, 1].set_ylabel("Frequency")

# subplot(1,1): capacity_factor histogram
axs[1, 1].hist(capacity_factor_dispatch, bins=40,
               range=[30, 100], color="blue")
axs[1, 1].grid(which="major", axis="both", linestyle="--")
axs[1, 1].set_xlabel("Capacity factor (%)")
axs[1, 1].set_ylabel("Frequency")

plt.tight_layout()
# plt.savefig("manuscript_figs/rts_pricetaker_full_year_4_subplots_scatter.pdf",
#                 format="pdf",
#                 bbox_inches="tight")
plt.savefig("manuscript_figs/rts_pricetaker_full_year_4_subplots_scatter.png",
        format="png", dpi=1000,
        bbox_inches="tight")

plt.show()
"""

# Generate histograms only for LMP and power produced
# Histogram for power produced
fig_1, ax_1 = plt.subplots()
ax_1.hist(capacity_factor_dispatch, bins=40,
          range=[30, 100], color="blue")
ax_1.grid(which="major", axis="both", linestyle="--")
ax_1.set_xlabel("Capacity factor (%)")
ax_1.set_ylabel("Hours in a year")
ax_1.text(60, 2000, "$P_{min} = 52.5$ MW", fontsize=12)
ax_1.text(60, 1900, "$P_{max} = 175$ MW", fontsize=12)
# fig_1.savefig("manuscript_figs/rts_pricetaker_full_year_power_hist.pdf",
#                 format="pdf",
#                 bbox_inches="tight")
fig_1.savefig("manuscript_figs/rts_price_taker_full_year_with_boil_eff_coal_50_dispatch.png",
        format="png", dpi=1000,
        bbox_inches="tight")

# Histogram for lmp
fig_2, ax_2 = plt.subplots()
ax_2.hist(price, bins=50, color="green")
ax_2.grid(which="major", axis="both", linestyle="--")
ax_2.set_xlabel("LMP ($/MWh)")
ax_2.set_ylabel("Hours in a year")
# fig_2.savefig("manuscript_figs/rts_pricetaker_full_year_price_hist.pdf",
#                 format="pdf",
#                 bbox_inches="tight")
fig_2.savefig("manuscript_figs/rts_price_taker_full_year_with_boil_eff_coal_50_lmp.png",
        format="png", dpi=1000,
        bbox_inches="tight")

plt.show()
