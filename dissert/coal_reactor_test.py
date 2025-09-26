# Verification test for biomass combustion stoichiometry

from pyomo.environ import ConcreteModel, SolverFactory
from idaes.models.unit_models import StoichiometricReactor
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog
from idaes.models.properties.modular_properties import GenericParameterBlock
from coal_comb_properties import configuration
from coal_comb_rp import MultiCombReactionParameterBlock

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.properties = GenericParameterBlock(**configuration)
m.fs.reaction = MultiCombReactionParameterBlock(property_package=m.fs.properties)


m.fs.react = StoichiometricReactor(
    property_package = m.fs.properties,
    reaction_package = m.fs.reaction,
    has_heat_of_reaction=True,
    has_heat_transfer=True, 
    has_pressure_change=False
)

M_coal = 100 # [g/s]
AFratio = 7.79
mw_air = 28.96 #[g/mol]
M_air = M_coal*AFratio
N_air = M_air/mw_air

mw_bm=configuration["components"]["coal"]["parameter_data"]["mw"][0]
N_bm = (M_coal/mw_bm)
N_total = N_bm + N_air


m.fs.react.config.reaction_package.rate_reaction_stoichiometry["Rcoal", "Sol", "uncombustible"].fix(0)
#mole_frac_comp spec
m.fs.react.inlet.mole_frac_comp[0,"O2"].fix(N_air*0.21/N_total)
m.fs.react.inlet.mole_frac_comp[0,"N2"].fix(N_air*0.79/N_total)
m.fs.react.inlet.mole_frac_comp[0,"CH4"].fix(1e-20)
m.fs.react.inlet.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.react.inlet.mole_frac_comp[0,"H2O"].fix(1e-20)
m.fs.react.inlet.mole_frac_comp[0,"coal"].fix(N_bm/N_total)
m.fs.react.inlet.mole_frac_comp[0,"uncombustible"].fix(1e-20)
m.fs.react.inlet.flow_mol.fix(N_total)
m.fs.react.inlet.temperature.fix(400)
m.fs.react.inlet.pressure.fix(101325)


m.fs.react.rate_reaction_extent[0,"RCH4"].fix(0)
m.fs.react.rate_reaction_extent[0,"Rcoal"].fix(N_bm)

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
