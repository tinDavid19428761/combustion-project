# Verification test for biomass combustion stoichiometry

from idaes.core.util.model_statistics import degrees_of_freedom

import idaes.logger as idaeslog

from pyomo.environ import ConcreteModel, SolverFactory
from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties import GenericParameterBlock
from bm_comb_properties import configuration
from bm_comb_rp import BMCombRxnParameterBlock
from idaes.models.unit_models import StoichiometricReactor

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

m.fs.properties = GenericParameterBlock(**configuration)
m.fs.reaction = BMCombRxnParameterBlock(property_package=m.fs.properties)
m.fs.react = StoichiometricReactor(
    property_package = m.fs.properties,
    reaction_package = m.fs.reaction,
    has_heat_of_reaction=True,
    has_heat_transfer=True, 
    has_pressure_change=False
)

#to embed into comprehensive combustion reactor unit model
# total_flow=1
# biomass_frac=0.1
M_bm = 100 # [g/s]
FAratio = 6.5
mw_air = 28.96 #[g/mol]
M_air = M_bm*FAratio
N_air = M_air/mw_air

ash_wt=0.02
w_bm = 0.09
mw_bm=configuration["components"]["biomass"]["parameter_data"]["mw"][0]
mw_ash=configuration["components"]["uncombustible"]["parameter_data"]["mw"][0]
N_bm = (M_bm/mw_bm)
N_ash = ash_wt*(1-w_bm)*M_bm/mw_ash
stoich_ash = N_ash/N_bm
N_total = N_bm + N_air

#adjusting biomass combustion stoichiometry for incombustible content:
new_co2 = 6*11/(11+stoich_ash)
new_h2o = 5*11/(11+stoich_ash)

m.fs.react.config.reaction_package.rate_reaction_stoichiometry["Rbiomass", "Sol", "uncombustible"].fix(stoich_ash)
m.fs.react.config.reaction_package.rate_reaction_stoichiometry["Rbiomass", "Vap", "CO2"].fix(new_co2)
m.fs.react.config.reaction_package.rate_reaction_stoichiometry["Rbiomass", "Vap", "H2O"].fix(new_h2o)

#mole_frac_comp spec
m.fs.react.inlet.mole_frac_comp[0,"O2"].fix(N_air*0.21/N_total)
m.fs.react.inlet.mole_frac_comp[0,"N2"].fix(N_air*0.79/N_total)
m.fs.react.inlet.mole_frac_comp[0,"CH4"].fix(1e-20)
m.fs.react.inlet.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.react.inlet.mole_frac_comp[0,"H2O"].fix(1e-20)
m.fs.react.inlet.mole_frac_comp[0,"biomass"].fix(N_bm/N_total)
m.fs.react.inlet.mole_frac_comp[0,"uncombustible"].fix(1e-20)
m.fs.react.inlet.flow_mol.fix(N_total)
m.fs.react.inlet.temperature.fix(400)
m.fs.react.inlet.pressure.fix(101325)


m.fs.react.rate_reaction_extent[0,"RCH4"].fix(0)
m.fs.react.rate_reaction_extent[0,"Rbiomass"].fix(N_bm)

# m.fs.react.heat_duty[0].fix(0)
m.fs.react.outlet.temperature.fix(300+273.15)

print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0
m.fs.react.initialize(outlvl=idaeslog.INFO)
solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)
m.fs.react.report()
# m.fs.react.display()

# print(value(m.fs.react.rate_reaction_extent))
