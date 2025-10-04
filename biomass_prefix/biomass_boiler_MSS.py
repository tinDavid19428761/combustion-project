""" 
Combustion Boiler Model with steam superheater. 
Modelled by adiabatic combustion reactor sending hot flue to boiler HX and superheater HX in counter-current to boiler water stream.

                                                            BlowDownWater
                                                                  ^
                                                                  |
                                                 +--------[Phase Separation]
                                                 |                |
                                                 |           BoilerWater
                                                 |                ^
                                                 |                |
Fuel+Air --->[Reactor]<|---> Hot Flue ---> [Superheater] ---> [Boiler] ---> Stack Flue
                       |                         |                ^
                       |                         \/               |
                       \/                    Superheated    BoilerFeedWater
                       Ash                      Steam             
                            
"""

__author__ = "David Dickson"

#Importing required pyomo and idaes components

import numpy as np
from pyomo.environ import (
    Constraint,
    Var,
    ConcreteModel,
    Expression,
    Param,
    Objective,
    SolverFactory,
    TransformationFactory,
    value,
    units as pyunits
)
from pyomo.network import Arc, SequentialDecomposition

#Todo add the four other unit operations
from idaes.models.unit_models import (
Mixer,
StoichiometricReactor,
Heater,
HeatExchanger,
Separator,
Flash
)
from idaes.models_extra.power_generation.unit_models.helm.phase_separator import HelmPhaseSeparator

from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
# Import idaes logger to set output levels
import idaes.logger as idaeslog
from  idaes.models.properties.modular_properties import GenericParameterBlock
from  biomass_comb_pp import configuration 
from  biomass_combustion_rp import BMCombReactionParameterBlock

#helmholtz import for water
from idaes.models.properties.general_helmholtz import (
        HelmholtzParameterBlock,
        HelmholtzThermoExpressions,
        AmountBasis,
        PhaseType,
        StateVars
    )
from idaes.models.unit_models.heat_exchanger import delta_temperature_amtd_callback, HX0DInitializer
from idaes.models.unit_models.separator import SplittingType, EnergySplittingType
from idaes.core.util.model_diagnostics import (
    DiagnosticsToolbox,
)
# import custom_combustion_reactor

import matplotlib.pyplot as plt

from custom_combustion_reactor import MultiCombReactor

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.biomass_properties = GenericParameterBlock(**configuration)
# m.fs.reaction_params = BMCombReactionParameterBlock(
#     property_package=m.fs.biomass_properties
# )

m.fs.steam_properties = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE,
        phase_presentation=PhaseType.LG,
        # state_vars=StateVars.TPX
    )

m.fs.fire_side = MultiCombReactor(
    property_package = m.fs.biomass_properties,
    # reaction_package = m.fs.reaction_params,
    has_heat_of_reaction=True,
    has_heat_transfer=True,
    has_pressure_change=False,
)

m.fs.superheater = HeatExchanger(
    delta_temperature_callback=delta_temperature_amtd_callback,
    hot_side_name="shell",
    cold_side_name="tube",
    shell={"property_package": m.fs.biomass_properties},
    tube={"property_package": m.fs.steam_properties}
)

m.fs.boiler_hx = HeatExchanger(
    delta_temperature_callback=delta_temperature_amtd_callback,
    hot_side_name="shell",
    cold_side_name="tube",
    shell={"property_package": m.fs.biomass_properties},
    tube={"property_package": m.fs.steam_properties}
)

m.fs.ash_sep = Separator(
    property_package = m.fs.biomass_properties,
    split_basis = SplittingType.phaseFlow,
    outlet_list = ["flue", "ash"]
)

m.fs.bdw_sep = HelmPhaseSeparator(
    property_package = m.fs.steam_properties,
)

#reactor flow sheet: feed-> reactor -> hot flue -> HX -> flue
m.fs.s01 = Arc(source=m.fs.fire_side.outlet,destination=m.fs.ash_sep.inlet)
m.fs.s02 = Arc(source=m.fs.ash_sep.flue,destination=m.fs.superheater.shell_inlet)
m.fs.s03 = Arc(source=m.fs.superheater.shell_outlet,destination=m.fs.boiler_hx.shell_inlet)
m.fs.s04 = Arc(source=m.fs.boiler_hx.tube_outlet,destination=m.fs.bdw_sep.inlet)
m.fs.s05 = Arc(source=m.fs.bdw_sep.vap_outlet,destination=m.fs.superheater.tube_inlet)
TransformationFactory("network.expand_arcs").apply_to(m)

#specifying reaction package parameters
# m.fs.reaction_params.h.fix(0.06) #h and w in dh_rxn calculation
# m.fs.reaction_params.w.fix(0.09)
# m.fs.reaction_params.rate_reaction_stoichiometry["R1","Sol","ash"]=0.02 #specify ash content

#reaction conversion constraint
# m.fs.fire_side.conversion = Var(initialize=1, bounds=(0,1))
# m.fs.fire_side.conversion_constraint = Constraint(
#     expr=m.fs.fire_side.conversion*m.fs.fire_side.inlet.mole_frac_comp[0,"biomass"]*m.fs.fire_side.inlet.flow_mol[0]
#     == (
#         m.fs.fire_side.inlet.flow_mol[0]*m.fs.fire_side.inlet.mole_frac_comp[0,"biomass"]
#         -m.fs.fire_side.outlet.flow_mol[0]*m.fs.fire_side.outlet.mole_frac_comp[0,"biomass"]
#     )
# )
m.fs.fire_side.conversion_R1.fix(1)
m.fs.fire_side.hcon.fix(0.07)
m.fs.fire_side.wcon.fix(0.1)
m.fs.fire_side.ohtc.fix(100)
m.fs.fire_side.surface_area.fix(0.1)
m.fs.fire_side.surface_temp.fix(60)
m.fs.fire_side.ash_mass_R1.fix(0.0)

#modelling Q_loss_casing_convection
# m.fs.fire_side.ohtc = Param(initialize=500, units=pyunits.J/pyunits.m**2/pyunits.K/pyunits.s, doc="boiler casing overall heat transfer coefficient")
# m.fs.fire_side.surface_area = Param(initialize=0.02, units=pyunits.m**2, doc="casing outer surface area")
# m.fs.fire_side.surface_temp = Param(initialize=55+273.12, units=pyunits.K, doc="outer skin temperature of boiler")

# m.fs.fire_side.heat_loss_constraint = Constraint(
#     expr=-m.fs.fire_side.heat_duty[0]==
#     (m.fs.fire_side.ohtc*m.fs.fire_side.surface_area*(m.fs.fire_side.outlet.temperature[0]-m.fs.fire_side.surface_temp))
# )

#reactor feed stream
m.fs.fire_side.inlet.mole_frac_comp[0,"N2"].fix(0.70)
m.fs.fire_side.inlet.mole_frac_comp[0,"O2"].fix(0.29)
m.fs.fire_side.inlet.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.fire_side.inlet.mole_frac_comp[0,"H2O"].fix(1e-20) 
m.fs.fire_side.inlet.mole_frac_comp[0,"biomass"].fix(0.01) 
m.fs.fire_side.inlet.mole_frac_comp[0,"ash"].fix(1e-20)
m.fs.fire_side.inlet.temperature.fix(300)
m.fs.fire_side.inlet.pressure.fix(101325)
m.fs.fire_side.inlet.flow_mol.fix(40)

#specifying ash separation
m.fs.ash_sep.split_fraction[0,"ash","Sol"].fix(1)
m.fs.ash_sep.split_fraction[0,"ash","Vap"].fix(0)

steam_press=101325

m.fs.boiler_hx.tube_inlet.enth_mol.fix(m.fs.steam_properties.htpx(p=steam_press*pyunits.Pa,T=300*pyunits.K))
m.fs.boiler_hx.tube_outlet.enth_mol.fix(m.fs.steam_properties.htpx(p=steam_press*pyunits.Pa,x=0.95)) # quality[x] = 1 - BlowDownRatio
m.fs.superheater.tube_outlet.enth_mol.fix(m.fs.steam_properties.htpx(p=steam_press*pyunits.Pa,T=400*pyunits.K))
m.fs.superheater.overall_heat_transfer_coefficient[0].fix(100)
m.fs.boiler_hx.tube_inlet.flow_mol.fix(20)
m.fs.boiler_hx.tube_inlet.pressure.fix(steam_press)
m.fs.boiler_hx.overall_heat_transfer_coefficient[0].fix(100)


#initialization routine
seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"
seq.options.tear_method = "Wegstein"
seq.options.iterLim = 3

G = seq.create_graph(m)
heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
order = seq.calculation_order(G)

#for identifying tear stream:
""" for o in heuristic_tear_set: #fs.s03
    print(o.name) """ 

tear_guesses = {
    "mole_frac_comp": {
        (0, "N2"): 0.67,
        (0, "O2"): 0.22,
        (0, "CO2"): 0.06,
        (0, "H2O"): 0.05,
        (0, "biomass"): 1e-20,
        (0, "ash"): 1e-20,
    },
    "flow_mol": {0: 40},
    "temperature": {0: 1000},
    "pressure": {0: 101325},
}

seq.set_guesses_for(m.fs.boiler_hx.shell_inlet, tear_guesses)

def function(unit):
    unit.initialize(outlvl=idaeslog.INFO)
    # print(degrees_of_freedom(unit))

print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0  
seq.run(m, function)

# m.fs.boiler_hx.overall_heat_transfer_coefficient[0].unfix()
# m.fs.boiler_hx.shell_outlet.temperature.fix(400)
m.fs.fire_side.ash_mass_R1.fix(0.01)
m.fs.fire_side.inlet.flow_mol.unfix()
solver=SolverFactory("ipopt")

m.fs.boiler_eff = Expression(
    expr = 100*(m.fs.superheater.heat_duty[0]+m.fs.boiler_hx.heat_duty[0])/(m.fs.fire_side.control_volume.rate_reaction_extent[0,"R1"]*-m.fs.fire_side.reaction_package.dh_rxn["R1"])
    )

m.fs.duty_to_steam = Expression(
    expr = (m.fs.superheater.heat_duty[0]+m.fs.boiler_hx.heat_duty[0])/1000/1000/1000 #GJ/s
)


# steady_states = [12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]
steady_states = list(range(2,30,2))
# steady_states = [20]
flue_temps = [340,370,400, 430]

efficiencies = np.zeros((len(flue_temps),len(steady_states)))
steam_duties = np.zeros((len(flue_temps),len(steady_states)))

""" Create MSS For-Loop Starting Here: """
for p,n in enumerate(flue_temps):
    m.fs.boiler_hx.shell_outlet.temperature.fix(n)
    for j,i in enumerate(steady_states):

        m.fs.boiler_hx.tube_inlet.flow_mol.fix(i)
        #pre-solve [actual] re-specification
        # m.fs.boiler_hx.shell_outlet.temperature.fix(450)
        # m.fs.boiler_hx.area.fix(20)

        # solver=SolverFactory("ipopt")
        print(degrees_of_freedom(m))
        assert degrees_of_freedom(m) == 0  
        status=solver.solve(m,tee=True)

        steam_duties[p][j] = value(m.fs.duty_to_steam)
        efficiencies[p][j] = value(m.fs.boiler_eff)
        print(f"{p}-{j}")

print(efficiencies)
print(steam_duties)


for i,j in enumerate(flue_temps):
    plt.plot(steam_duties[i],efficiencies[i],label=f"{(j-273.15):.2f} C ")
plt.ylim(80,100)
plt.xlabel('Steam Generation Duty [GJ/s]')
plt.ylabel('Boiler System Efficiency [%]')
plt.legend()
plt.show()


#results
m.fs.fire_side.report()
m.fs.superheater.report()
m.fs.boiler_hx.report()

print(f"    Boiler Efficiency: {value(m.fs.boiler_eff):.2f}%")
# print(f"    casing heat loss:{value(m.fs.fire_side.heat_duty[0]):.2f}J/s")

