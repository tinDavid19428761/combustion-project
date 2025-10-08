"""
custom combustion reactor unit model test:::
"""
#Importing required pyomo and idaes components
from pyomo.environ import (
    Constraint,
    Var,
    ConcreteModel,
    Expression,
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
)

from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
# Import idaes logger to set output levels
import idaes.logger as idaeslog
from idaes.models.properties.modular_properties import GenericParameterBlock
from  biomass_comb_pp import configuration 

#helmholtz import for water
from idaes.models.properties.general_helmholtz import (
        HelmholtzParameterBlock,
        HelmholtzThermoExpressions,
        AmountBasis,
        PhaseType,
    )

# from custom_combustion_airspec import MultiCombReactor
from custom_combustion_reactor import MultiCombReactor

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.biomass_properties = GenericParameterBlock(**configuration)

m.fs.steam_properties = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE,
        phase_presentation=PhaseType.LG,
    )

from idaes.core.util.model_statistics import degrees_of_freedom



m.fs.R101 = MultiCombReactor(
    property_package = m.fs.biomass_properties,
    has_heat_of_reaction=True,
    has_heat_transfer=True,
    has_pressure_change=False,
)




#air_spec fix var constraints:
# m.fs.R101.excess_air_percent_Rcoal.fix(50)
# m.fs.R101.air_fuel_mass_ratio_R1.fix(10)
# m.fs.R101.mass_flow_kg_coal.fix(1.5)
# m.fs.R101.mass_flow_kg_biomass.fix(1)



m.fs.R101.conversion_Rcoal.fix(1)
m.fs.R101.ash_mass_Rcoal.fix(0.0)



m.fs.R101.conversion_R1.fix(1)
m.fs.R101.hcon.fix(0.06)
m.fs.R101.wcon.fix(0.09)
m.fs.R101.ohtc.fix(100)
m.fs.R101.surface_temp.fix(60)
m.fs.R101.ash_mass_R1.fix(0.0)


m.fs.R101.inlet.mole_frac_comp[0,"N2"].fix(0.7)
m.fs.R101.inlet.mole_frac_comp[0,"O2"].fix(0.3-0.1-0.01)
m.fs.R101.inlet.mole_frac_comp[0,"biomass"].fix(0.01) 
m.fs.R101.inlet.mole_frac_comp[0,"coal"].fix(0.1) 

m.fs.R101.inlet.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.R101.inlet.mole_frac_comp[0,"H2O"].fix(1e-20) 


m.fs.R101.inlet.mole_frac_comp[0,"ash"].fix(1e-20) 
m.fs.R101.inlet.temperature.fix(300)
m.fs.R101.inlet.pressure.fix(101325)

#remove when using airspec
m.fs.R101.inlet.flow_mol.fix(10)

m.fs.R101.outlet.temperature.fix(400)

print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0

m.fs.R101.initialize(outlvl=idaeslog.INFO)

m.fs.R101.ash_mass_R1.fix(0.03)
m.fs.R101.ash_mass_Rcoal.fix(0.0)


solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)

m.fs.R101.report()
print(value(m.fs.R101.reaction_package.rate_reaction_stoichiometry["R1", "Sol", "biomass"]))
print(value(m.fs.R101.reaction_package.rate_reaction_stoichiometry["R1", "Sol", "ash"]))
