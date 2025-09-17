""" 
Combustion Boiler Model with steam superheater. 
Modelled by adiabatic combustion reactor sending hot flue to boiler HX and superheater HX in counter-current to boiler water stream.

                                                            BlowDownWater
                                                                  ^
                                                                  |
                                                 +--------[Phase Separation]
                                                 |                |
                                                 |           BoilerWater   BoilerFeedWater
                                                 |                ^             |
                                                 |                |             V
Fuel+Air --->[Reactor]<|---> Hot Flue ---> [Superheater] ---> [Boiler] ---> [Economizer] ---> Stack Flue
                       |                         |                ^             |
                       |                         V                |             |
                       V                    Superheated           |--heated feed--
                      Ash                      Steam              
                            
"""

__author__ = "David Dickson"

#Importing required pyomo and idaes components
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
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent, check_units_equivalent
#Todo add the four other unit operations
from idaes.models.unit_models import (
Mixer,
# StoichiometricReactor,
Heater,
HeatExchanger,
Separator,
Flash
)
# from custom_stoichiometric_reactor import StoichiometricReactor
from custom_combustion_reactor import MultiCombReactor

from idaes.models_extra.power_generation.unit_models.helm.phase_separator import HelmPhaseSeparator
from idaes.models_extra.power_generation.unit_models.boiler_heat_exchanger import (
    BoilerHeatExchanger,
    TubeArrangement,
    HeatExchangerFlowPattern,
)

from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
# Import idaes logger to set output levels
import idaes.logger as idaeslog
from idaes.models.properties.modular_properties import GenericParameterBlock
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
# from idaes.models_extra.power_generation.properties.flue_gas_ideal import FlueGasParameterBlock
from random.parameter_block_flue_gas_ideal import FlueGasParameterBlock


m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.biomass_properties = GenericParameterBlock(**configuration)

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

# m.fs.superheater = HeatExchanger(
#     delta_temperature_callback=delta_temperature_amtd_callback,
#     hot_side_name="shell",
#     cold_side_name="tube",
#     shell={"property_package": m.fs.biomass_properties,
#            "has_pressure_change": False},
#     tube={"property_package": m.fs.steam_properties,
#           "has_pressure_change": False}
# )

m.fs.boiler_hx = HeatExchanger(
    delta_temperature_callback=delta_temperature_amtd_callback,
    hot_side_name="shell",
    cold_side_name="tube",
    shell={"property_package": m.fs.biomass_properties},
    tube={"property_package": m.fs.steam_properties}
)

m.fs.econ = HeatExchanger(
    delta_temperature_callback=delta_temperature_amtd_callback,
    hot_side_name="shell",
    cold_side_name="tube",
    shell={"property_package": m.fs.biomass_properties},
    tube={"property_package": m.fs.steam_properties}
)

m.fs.ash_sep = Separator(
    property_package = m.fs.biomass_properties,
    split_basis = SplittingType.phaseFlow,
    outlet_list = ["flue", "ash_out"]
)

m.fs.bdw_sep = HelmPhaseSeparator(
    property_package = m.fs.steam_properties,
)

#reactor flow sheet: feed-> reactor -> hot flue -> HX -> flue
m.fs.s01 = Arc(source=m.fs.fire_side.outlet,destination=m.fs.ash_sep.inlet)
m.fs.s02 = Arc(source=m.fs.ash_sep.flue,destination=m.fs.boiler_hx.shell_inlet)
# m.fs.s03 = Arc(source=m.fs.superheater.shell_outlet,destination=m.fs.boiler_hx.shell_inlet)
m.fs.s04 = Arc(source=m.fs.boiler_hx.tube_outlet,destination=m.fs.bdw_sep.inlet)
# m.fs.s05 = Arc(source=m.fs.bdw_sep.vap_outlet,destination=m.fs.superheater.tube_inlet)
m.fs.s06 = Arc(source=m.fs.econ.tube_outlet, destination=m.fs.boiler_hx.tube_inlet)
m.fs.s07 = Arc(source=m.fs.boiler_hx.shell_outlet, destination=m.fs.econ.shell_inlet)

TransformationFactory("network.expand_arcs").apply_to(m)

#specifying reaction package variables
m.fs.fire_side.reaction_package.h.fix(0.06) #h and w in dh_rxn calculation
m.fs.fire_side.reaction_package.w.fix(0.09)
m.fs.fire_side.reaction_package.rate_reaction_stoichiometry["Rbiomass","Sol","ash"].unfix()
m.fs.fire_side.ash_mass.fix(0.03)

#specifying combustion reactor variables
m.fs.fire_side.conversion["Rbiomass"].fix(1)
m.fs.fire_side.conversion["RCH4"].fix(1)
m.fs.fire_side.heat_duty[0].fix(-1000)
m.fs.fire_side.surface_area.fix(0.1)
m.fs.fire_side.surface_temp.fix(55+273.15)

#specifying feed stream variables
m.fs.fire_side.inlet.mole_frac_comp[0,"N2"].fix(0.5)
m.fs.fire_side.inlet.mole_frac_comp[0,"O2"].fix(0.39)
m.fs.fire_side.inlet.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.fire_side.inlet.mole_frac_comp[0,"H2O"].fix(1e-20) 
m.fs.fire_side.inlet.mole_frac_comp[0,"CO"].fix(1e-20) 
m.fs.fire_side.inlet.mole_frac_comp[0,"biomass"].fix(0.01) 
m.fs.fire_side.inlet.mole_frac_comp[0,"ash"].fix(1e-20)
m.fs.fire_side.inlet.mole_frac_comp[0,"CH4"].fix(0.1)
m.fs.fire_side.inlet.temperature.fix(300)
m.fs.fire_side.inlet.pressure.fix(101325)
m.fs.fire_side.inlet.flow_mol.fix(40)

#specifying ash separator variables
m.fs.ash_sep.split_fraction[0,"ash_out","Sol"].fix(1)
m.fs.ash_sep.split_fraction[0,"ash_out","Vap"].fix(0)

#specifying water/steam stream variables

m.fs.boiler_hx.tube_outlet.enth_mol.fix(m.fs.steam_properties.htpx(p=101325*pyunits.Pa,x=0.95))
m.fs.boiler_hx.overall_heat_transfer_coefficient[0].fix(1000)

#specifying economizer
m.fs.econ.tube_inlet.flow_mol.fix(20)
m.fs.econ.tube_inlet.pressure.fix(101325)
m.fs.econ.tube_inlet.enth_mol.fix(m.fs.steam_properties.htpx(p=101325*pyunits.Pa,T=300*pyunits.K))
# m.fs.boiler_hx.tube_inlet.enth_mol.fix(m.fs.steam_properties.htpx(p=101325*pyunits.Pa,T=310*pyunits.K))           
m.fs.econ.shell_outlet.temperature.fix(400)
m.fs.econ.overall_heat_transfer_coefficient[0].fix(1000)
# m.fs.econ.area.fix(20)

# m.fs.boiler_hx.shell_outlet.temperature.fix(420)


#initialization routine
seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"
seq.options.tear_method = "Wegstein"
seq.options.iterLim = 3

G = seq.create_graph(m)
heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
order = seq.calculation_order(G)

# for identifying tear stream:
for o in heuristic_tear_set: #fs.s03, fs.s07
    print(o.name)

tear_guesses = {
    "pre_boiler":{
    "mole_frac_comp": {
        (0, "N2"): 0.67,
        (0, "O2"): 0.22,
        (0, "CO2"): 0.06,
        (0, "H2O"): 0.05,
        (0, "CO"): 1e-20,
        (0, "CH4"): 1e-20,
        (0, "biomass"): 1e-20,
        (0, "ash"): 1e-20,
    },
    "flow_mol": {0: 40},
    "temperature": {0: 900},
    "pressure": {0: 101325},
    },
    "post_boiler":{
    "mole_frac_comp": {
        (0, "N2"): 0.67,
        (0, "O2"): 0.22,
        (0, "CO2"): 0.06,
        (0, "H2O"): 0.05,
        (0, "CO"): 1e-20,
        (0, "CH4"): 1e-20,
        (0, "biomass"): 1e-20,
        (0, "ash"): 1e-20,
    },
    "flow_mol": {0: 40},
    "temperature": {0: 420},
    "pressure": {0: 101325},
    },

}

# seq.set_guesses_for(m.fs.boiler_hx.shell_inlet, tear_guesses["pre_boiler"])
seq.set_guesses_for(m.fs.boiler_hx.shell_outlet, tear_guesses["post_boiler"])

def function(unit):
    print(degrees_of_freedom(unit))
    unit.initialize(outlvl=idaeslog.INFO)
    unit.report()
print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0

seq.run(m, function)

# dt = DiagnosticsToolbox(m)
# dt.report_structural_issues()
# dt.display_underconstrained_set()
# dt.display_components_with_inconsistent_units()

#pre-solve [actual] re-specification
m.fs.fire_side.inlet.flow_mol.unfix()

solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)

m.fs.boiler_eff = Expression(
    expr = (m.fs.superheater.heat_duty[0]+m.fs.boiler_hx.heat_duty[0])/(m.fs.fire_side.control_volume.rate_reaction_extent[0,"Rbiomass"]*-m.fs.fire_side.reaction_package.dh_rxn["Rbiomass"])
)
m.fs.biomass_mass_flow = Expression(#units g/s 
    expr = (m.fs.fire_side.inlet.mole_frac_comp[0,"biomass"])*(m.fs.fire_side.inlet.flow_mol[0])*m.fs.biomass_properties.biomass.mw
)

#results
m.fs.fire_side.report()
m.fs.superheater.report()
m.fs.boiler_hx.report()
m.fs.bdw_sep.report()
m.fs.ash_sep.report()
print(f"    Boiler Efficiency: {value(m.fs.boiler_eff)*100:.2f}%")
print(f"    casing heat loss:{value(m.fs.fire_side.heat_duty[0]):.2f} J/s")
print(f"    biomass demand: {value(m.fs.biomass_mass_flow):.3f} g/s")
print(value(m.fs.fire_side.reaction_package.rate_reaction_stoichiometry["Rbiomass","Sol","ash"]))


