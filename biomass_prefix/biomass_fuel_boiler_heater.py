# -*- coding: utf-8 -*-
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
from  biomass_combustion_rp import BMCombReactionParameterBlock

#helmholtz import for water
from idaes.models.properties.general_helmholtz import (
        HelmholtzParameterBlock,
        HelmholtzThermoExpressions,
        AmountBasis,
        PhaseType,
    )

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
    )

from idaes.core.util.model_statistics import degrees_of_freedom



m.fs.R101 = MultiCombReactor(
    property_package = m.fs.biomass_properties,
    # reaction_package = m.fs.reaction_params,
    has_heat_of_reaction=True,
    has_heat_transfer=True,
    has_pressure_change=False,
)



flowTotal = 100
bm_frac=0.01

m.fs.R101.conversion_Rcoal.fix(1)
m.fs.R101.ash_mass_Rcoal.fix(0.0)
m.fs.R101.dh_rxn_Rcoal.fix(-284675.1254)


#appended variables by custom_combustion_reactor
m.fs.R101.conversion_R1.fix(1)
m.fs.R101.hcon.fix(0.06)
m.fs.R101.wcon.fix(0.09)
m.fs.R101.ohtc.fix(100)
# m.fs.R101.surface_area.fix(1)
m.fs.R101.surface_temp.fix(60)
m.fs.R101.ash_mass_R1.fix(0.0)
# m.fs.R101.heat_duty.fix(-100000)


m.fs.R101.inlet.mole_frac_comp[0,"N2"].fix(0.7)
m.fs.R101.inlet.mole_frac_comp[0,"O2"].fix(0.28)
m.fs.R101.inlet.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.R101.inlet.mole_frac_comp[0,"H2O"].fix(1e-20) 
m.fs.R101.inlet.mole_frac_comp[0,"biomass"].fix(bm_frac) 
m.fs.R101.inlet.mole_frac_comp[0,"coal"].fix(bm_frac) 
m.fs.R101.inlet.mole_frac_comp[0,"ash"].fix(1e-20) 
m.fs.R101.inlet.temperature.fix(300)
m.fs.R101.inlet.pressure.fix(101325)
m.fs.R101.inlet.flow_mol.fix(flowTotal)


m.fs.R101.outlet.temperature.fix(400)
# m.fs.R101.rate_reaction_extent[0,"R1"].fix(extentR1)


print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0
m.fs.R101.initialize(outlvl=idaeslog.INFO)

m.fs.R101.outlet.temperature.unfix()
m.fs.R101.surface_area.fix(0.1)
m.fs.R101.ash_mass_R1.fix(0.02)
m.fs.R101.ash_mass_Rcoal.fix(0.03)


solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)

m.fs.R101.report()
print(value(m.fs.R101.reaction_package.rate_reaction_stoichiometry["R1", "Sol", "biomass"]))
print(value(m.fs.R101.reaction_package.rate_reaction_stoichiometry["R1", "Sol", "ash"]))
# print(f"{value(m.fs.R101.reaction_package.rate_reaction_stoichiometry[("R1", "Vap", "H2O")]):.2f}")
# print(f"{value(m.fs.R101.reaction_package.rate_reaction_stoichiometry["R1", "Sol", "ash"]):.8f}")
# print(f"{value(m.fs.R101.dh_rxn_R1):.2f}")

# print(f"{value(m.fs.R101.reaction_package.rate_reaction_stoichiometry["R1", "Vap", "N2"]):.2f}")
# print(f"{value(m.fs.R101.reaction_package.rate_reaction_stoichiometry["R1", "Vap", "O2"]):.2f}")


# m.fs.H101 = Heater(
#     property_package = m.fs.steam_properties,
#     has_pressure_change=False,
#     has_phase_equilibrium=False,
# )

# #steam generation specification
# steam_pressure = 101325
# T_water_in=300
# T_water_out=380
# enth_in=m.fs.steam_properties.htpx(p=steam_pressure * pyunits.Pa, T= T_water_in * pyunits.K)
# enth_out=m.fs.steam_properties.htpx(p=steam_pressure * pyunits.Pa, T= T_water_out * pyunits.K)

# #steam_gen for initialization purposes
# initSteam = 10

# m.fs.H101.inlet.flow_mol.fix(initSteam)
# m.fs.H101.inlet.enth_mol.fix(enth_in)
# m.fs.H101.inlet.pressure.fix(steam_pressure)
# m.fs.H101.heat_duty.fix(12345)

# m.fs.H101.initialize(outlvl=idaeslog.INFO)

# #re-specifying .fix's
# m.fs.H101.inlet.flow_mol.unfix()
# m.fs.H101.outlet.enth_mol.fix(enth_out)

# m.fs.H101.heat_duty.fix(-value(m.fs.R101.heat_duty[0]))

# #solve fow steam generation flow
# solver=SolverFactory("ipopt")
# status=solver.solve(m,tee=True)
# # print(value(m.fs.R101.heat_duty[0]))
# print(value(m.fs.H101.heat_duty[0]))

# m.fs.boiler_eff = Expression(
#         expr = 100*(m.fs.H101.heat_duty[0])/(m.fs.R101.rate_reaction_extent[0,"R1"]*-m.fs.reaction_params.dh_rxn["R1"])
#     )

# m.fs.R101.report()
# m.fs.H101.report()
# print(value(m.fs.boiler_eff))
