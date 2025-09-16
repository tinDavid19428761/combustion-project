""" 
Combustion Boiler Model with steam superheater. 
Modelled by adiabatic combustion reactor sending hot flue to boiler HX and superheater HX in counter-current to boiler water stream.




                                       BoilerWater                 
                                           ^            
                                           |  
Fuel+Air --->[Reactor]<|---> Hot Flue -[Boiler] ---> Stack Flue
                       |                   ^             
                       |                   |             
                       V             BoilerFeedWater
                      Ash(omit)                        
                            
"""

#stoich_reactor_test
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
# from pytest import approx

#Todo add the four other unit operations
from idaes.models.unit_models import (
Mixer,
# StoichiometricReactor,
Heater,
HeatExchanger,
Separator,
Flash,
Heater
)
from idaes.models_extra.power_generation.unit_models.helm.phase_separator import HelmPhaseSeparator

from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
# Import idaes logger to set output levels
import idaes.logger as idaeslog
from idaes.models.properties.modular_properties import GenericParameterBlock
# from  biomass_comb_pp import configuration 
from black_liquor_pp import configuration
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
from bl_combustion_reactor import MultiCombReactor #has rxn pkg included

# from custom_combustion_reactor import MultiCombReactor
import unittest

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.bl_properties = GenericParameterBlock(**configuration)



m.fs.steam_properties = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE,
        phase_presentation=PhaseType.LG,
    )

m.fs.R101 = MultiCombReactor(
    property_package = m.fs.bl_properties,
    # reaction_package = m.fs.reaction_params,
    has_heat_of_reaction=True,
    has_heat_transfer=True, 
    has_pressure_change=False
)

m.fs.H101 = Heater(
    property_package = m.fs.steam_properties,
)

# m.fs.mix = Mixer(
#     property_package = m.fs.bl_properties,
#     inlet_list = ["air","fuel"]
# )

# m.fs.s01 = Arc(source=m.fs.mix.outlet, destination=m.fs.R101.inlet)
# TransformationFactory("network.expand_arcs").apply_to(m)


#case study inputs
#flows
f = {
    # "air_flow": 976,
    "air_flow": 1670, #finding: will not feasible converge given limited oxygen. solution: increase air (strays from case study) or decrease stoich air in rxnpkg
    "bl_flow": 998,
    "gas_flow": 0.1717
}

flow_total = f["bl_flow"]+f["gas_flow"]+f["air_flow"]
ar = f["air_flow"]/(flow_total)

m.fs.R101.inlet.mole_frac_comp[0,"N2"].fix(0.79*ar)
m.fs.R101.inlet.mole_frac_comp[0,"O2"].fix(0.21*ar)
m.fs.R101.inlet.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.R101.inlet.mole_frac_comp[0,"H2O"].fix(1e-20)
m.fs.R101.inlet.mole_frac_comp[0,"BL"].fix(f["bl_flow"]/flow_total)
m.fs.R101.inlet.mole_frac_comp[0,"uncombustible"].fix(1e-20)
m.fs.R101.inlet.mole_frac_comp[0,"CH4"].fix(f["gas_flow"]/flow_total)
m.fs.R101.inlet.temperature.fix(100+273.15) #guess, switch to mixer for more accurate
m.fs.R101.inlet.pressure.fix(101325) #
m.fs.R101.inlet.flow_mol.fix(flow_total/10)

# #air stream
# m.fs.mix.air.mole_frac_comp[0,"N2"].fix(0.79)
# m.fs.mix.air.mole_frac_comp[0,"O2"].fix(0.21)
# m.fs.mix.air.mole_frac_comp[0,"CO2"].fix(1e-20)
# m.fs.mix.air.mole_frac_comp[0,"H2O"].fix(1e-20)
# m.fs.mix.air.mole_frac_comp[0,"BL"].fix(1e-20)
# m.fs.mix.air.mole_frac_comp[0,"uncombustible"].fix(1e-20)
# m.fs.mix.air.mole_frac_comp[0,"CH4"].fix(1e-20)
# m.fs.mix.air.temperature.fix(33+273.15)
# m.fs.mix.air.pressure.fix(101325)
# m.fs.mix.air.flow_mol.fix(f["air_flow"]/10) # /10 for keeping in initialisation bounds

# #fuel stream
# m.fs.mix.fuel.mole_frac_comp[0,"N2"].fix(1e-20)
# m.fs.mix.fuel.mole_frac_comp[0,"O2"].fix(1e-20)
# m.fs.mix.fuel.mole_frac_comp[0,"CO2"].fix(1e-20)
# m.fs.mix.fuel.mole_frac_comp[0,"H2O"].fix(1e-20)
# # m.fs.mix.fuel.mole_frac_comp[0,"BL"].fix(f["bl_flow"]/fuel_total)
# m.fs.mix.fuel.mole_frac_comp[0,"BL"].fix(1)
# m.fs.mix.fuel.mole_frac_comp[0,"uncombustible"].fix(1e-20)
# # m.fs.mix.fuel.mole_frac_comp[0,"CH4"].fix(f["gas_flow"]/fuel_total)
# m.fs.mix.fuel.mole_frac_comp[0,"CH4"].fix(1e-20)
# m.fs.mix.fuel.temperature.fix(123.6+273.15)
# m.fs.mix.fuel.pressure.fix(101325)
# m.fs.mix.fuel.flow_mol.fix(fuel_total/10)



m.fs.R101.conversion_Rbl.fix(1)
m.fs.R101.conversion_RCH4.fix(0)

m.fs.R101.dh_rxn_RCH4.fix(-802125) #must fix with actual numbers not just .fix()
m.fs.R101.dh_rxn_Rbl.fix(-135150)

m.fs.R101.outlet.temperature.fix(250+273.15)




m.fs.R101.initialize(outlvl=idaeslog.INFO)


m.fs.H101.inlet.flow_mol.fix(300)
m.fs.H101.inlet.enth_mol.fix(m.fs.steam_properties.htpx(p=101325*pyunits.Pa,T=300*pyunits.K))
m.fs.H101.inlet.pressure.fix(101325)

m.fs.H101.outlet.enth_mol.fix(m.fs.steam_properties.htpx(p=101325*pyunits.Pa,T=(400+273.15)*pyunits.K))




# m.fs.mix.air.flow_mol.fix(f["air_flow"]) #return to full flow values for solve
# m.fs.mix.fuel.flow_mol.fix(fuel_total)
print(degrees_of_freedom(m.fs.H101))
m.fs.H101.initialize(outlvl=idaeslog.INFO)

# m.fs.H101.inlet.flow_mol.fix(217.36)
m.fs.H101.inlet.flow_mol.unfix()
m.fs.H101.heat_duty.fix(-value(m.fs.R101.heat_duty[0]))
# m.fs.R101.heat_duty.fix(-value(m.fs.H101.heat_duty[0]))
# m.fs.R101.outlet.temperature.unfix() #solving for flue temp




print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0
solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)
m.fs.R101.report()
m.fs.H101.report()
# m.fs.R101.display()

