##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Simple rankine cycle model. Has couple of options:
1. Recover waste heat after turbine to mimic feed water heater integration
2. Option to include boiler efficiency which is a linear fit f(capacity factor)

if no heat recovery, the flowsheet is as follows:
    Boiler --> Turbine --> Condenser --> Pump --> Boiler

if heat_recovery, the flowsheet is as follows:
    Boiler --> Turbine --> pre-condenser(- Q_recovered) --> Condenser -->
    Pump --> Feed water heater(+ Q_recovered) --> Boiler

Note:
* Boiler and condenser are simple heater blocks
* IAPWS95 for water and steam properties
"""

__author__ = "Jaffer Ghouse"


# Import Pyomo libraries
from pyomo.environ import ConcreteModel, units, Var, \
    TransformationFactory, value, Block, Expression, Constraint, Param, \
    Objective
from pyomo.network import Arc

from matplotlib import pyplot as plt
from idaes.core.util import to_json, from_json

from simple_rankine_cycle import create_model, set_inputs, initialize_model, \
    close_flowsheet_loop, add_capital_cost, add_operating_cost


def steady_state_price_taker(heat_recovery=False,
                             calc_boiler_eff=False,
                             include_shutdown=False,
                             p_max_lower_bound=10,
                             p_max_upper_bound=300,
                             p_min_multiplier=0.3,
                             capital_payment_years=5,
                             plant_lifetime=20,
                             power_demand=None, lmp=None,
                             lmp_weights=None,
                             coal_price=None):
    """This method sets up the stochastic optimization problem that sets up a
    steady-state, pricetaker problem.

    Returns:
        [type]: [description]
    """
    m = ConcreteModel()

    # Create capex plant
    m.cap_fs = create_model(
        heat_recovery=heat_recovery,
        capital_fs=True, calc_boiler_eff=False)
    m.cap_fs = set_inputs(m.cap_fs)
    m.cap_fs = initialize_model(m.cap_fs)
    m.cap_fs = close_flowsheet_loop(m.cap_fs)
    m.cap_fs = add_capital_cost(m.cap_fs)

    # capital cost (M$/yr)
    cap_expr = m.cap_fs.fs.capital_cost*1e6/capital_payment_years

    # Create opex plant
    op_expr = 0
    rev_expr = 0

    for i in range(len(lmp)):

        print()
        print("Creating instance ", i)
        if not calc_boiler_eff:
            op_fs = create_model(
                heat_recovery=heat_recovery,
                capital_fs=False,
                calc_boiler_eff=False)

            # Set model inputs for the capex and opex plant
            op_fs = set_inputs(op_fs)

            if i == 0:
                # Initialize the capex and opex plant
                op_fs = initialize_model(op_fs)

                # save model state after initializing the first instance
                to_json(op_fs.fs, fname="initialized_state.json.gz",
                        gz=True, human_read=True)
            else:
                # Initialize the capex and opex plant
                from_json(op_fs.fs, fname="initialized_state.json.gz", gz=True)

            # Closing the loop in the flowsheet
            op_fs = close_flowsheet_loop(op_fs)

            op_fs = add_operating_cost(op_fs, coal_price=coal_price)

            op_expr += lmp_weights[i]*op_fs.fs.operating_cost
            rev_expr += lmp_weights[i]*lmp[i]* \
                op_fs.fs.net_cycle_power_output*1e-6

            # Add inequality constraint linking net power to cap_ex
            # operating P_min <= 30% of design P_max
            op_fs.fs.eq_min_power = Constraint(
                expr=op_fs.fs.net_cycle_power_output >=
                p_min_multiplier*m.cap_fs.fs.net_cycle_power_output)
            # operating P_max = design P_max
            op_fs.fs.eq_max_power = Constraint(
                expr=op_fs.fs.net_cycle_power_output <=
                m.cap_fs.fs.net_cycle_power_output)
        else:
            op_fs = create_model(
                heat_recovery=heat_recovery,
                capital_fs=False,
                calc_boiler_eff=True)

            # Set model inputs for the capex and opex plant
            op_fs = set_inputs(op_fs)

            # Fix the p_max of op_fs to p of cap_fs for initialization
            op_fs.fs.net_power_max.fix(
                value(m.cap_fs.fs.net_cycle_power_output))

            if i == 0:
                # Initialize the capex and opex plant
                op_fs = initialize_model(op_fs)

                # save model state after initializing the first instance
                to_json(op_fs.fs, fname="initialized_state.json.gz",
                        gz=True, human_read=True)
            else:
                # Initialize the capex and opex plant
                from_json(op_fs.fs, fname="initialized_state.json.gz", gz=True)

            # Closing the loop in the flowsheet
            op_fs = close_flowsheet_loop(op_fs)
            op_fs = add_operating_cost(op_fs, coal_price=coal_price)

            op_expr += lmp_weights[i]*op_fs.fs.on_off*op_fs.fs.operating_cost
            rev_expr += lmp_weights[i]*float(lmp[i])*op_fs.fs.on_off*op_fs.\
                fs.net_cycle_power_output*1e-6

            # Unfix op_fs p_max and set constraint linking that to cap_fs p_max
            op_fs.fs.net_power_max.unfix()
            op_fs.fs.eq_p_max = Constraint(
                expr=op_fs.fs.net_power_max ==
                m.cap_fs.fs.net_cycle_power_output*1e-6
            )

            # Add inequality constraint linking net power to cap_ex
            # operating P_min <= 30% of design P_max
            op_fs.fs.eq_min_power = Constraint(
                expr=op_fs.fs.net_cycle_power_output >=
                p_min_multiplier*m.cap_fs.fs.net_cycle_power_output)
            # operating P_max = design P_max
            op_fs.fs.eq_max_power = Constraint(
                expr=op_fs.fs.net_cycle_power_output <=
                m.cap_fs.fs.net_cycle_power_output)

            # Unfix the on-off variable
            if include_shutdown:
                op_fs.fs.on_off.unfix()
                op_fs.fs.on_off.setlb(0)
                op_fs.fs.on_off.setub(1)


        # only if power demand is given
        if power_demand is not None:
            op_fs.fs.eq_max_produced = Constraint(
                expr=op_fs.fs.net_cycle_power_output <=
                power_demand[i]*1e6)

        op_fs.fs.boiler.inlet.flow_mol[0].unfix()

        # Set bounds for the flow
        op_fs.fs.boiler.inlet.flow_mol[0].setlb(1)
        # op_fs.fs.boiler.inlet.flow_mol[0].setub(25000)

        setattr(m, 'scenario_{}'.format(i), op_fs)

    # Expression for total cap and op cost - $
    m.total_cost = Expression(
        expr=plant_lifetime*op_expr + capital_payment_years*cap_expr)

    # Expression for total revenue
    m.total_revenue = Expression(
        expr=plant_lifetime*rev_expr)

    # Objective $
    m.obj = Objective(
        expr=-(m.total_revenue - m.total_cost))

    # Unfixing the boiler inlet flowrate for capex plant
    m.cap_fs.fs.boiler.inlet.flow_mol[0].unfix()

    # Setting bounds for the capex plant flowrate
    m.cap_fs.fs.boiler.inlet.flow_mol[0].setlb(5)

    # Setting bounds for net cycle power output for the capex plant
    m.cap_fs.fs.eq_min_power = Constraint(
        expr=m.cap_fs.fs.net_cycle_power_output >= p_max_lower_bound*1e6)

    m.cap_fs.fs.eq_max_power = Constraint(
        expr=m.cap_fs.fs.net_cycle_power_output <=
        p_max_upper_bound*1e6)

    return m

