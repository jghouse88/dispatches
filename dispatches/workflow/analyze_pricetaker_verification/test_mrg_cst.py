import pyomo.environ as pyomo
from pyomo.common.fileutils import this_file_dir
import sys, os, json
import pandas as pd
import numpy as np

# add rankine cycle code to path
sys.path.append(os.path.join(this_file_dir(),"../../models/simple_case"))
import simple_rankine_cycle as src


pmax = 177.5
output_value1 = 177.5
coal_price=30.0

# operating cost at pmax
m1 = src.square_problem(
    heat_recovery=True,
    capital_fs=False,
    calc_boiler_eff=True,
    p_max=pmax,     # The pmax used for efficiency calculations
    net_power=output_value1, # The operating output (used as nameplate for capex)
    coal_price=coal_price
)
op_cost1 = pyomo.value(m1.fs.operating_cost)/output_value1

# operating cost at pmin=0.15*pmax
output_value2 = 177.5*0.15
m2 = src.square_problem(
    heat_recovery=True,
    capital_fs=False,
    calc_boiler_eff=True,
    p_max=pmax,     # The pmax used for efficiency calculations
    net_power=output_value2, # The operating output (used as nameplate for capex)
    coal_price=coal_price
)
op_cost2 = pyomo.value(m2.fs.operating_cost)/output_value2

# take average of two operating costs to get the average marginal cost
mrg_cost = 0.5*(op_cost1 + op_cost2)