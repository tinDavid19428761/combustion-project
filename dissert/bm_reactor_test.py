#modular_bm_reactor_test

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
StoichiometricReactor,
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
from idaes.models.properties.modular_properties import GenericParameterBlock, GenericReactionParameterBlock
# from  biomass_comb_pp import configuration 


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

# from custom_combustion_reactor import MultiCombReactor
import unittest

from biomass_comb_properties import configuration
from bm_comb_rp import MultiCombReactionParameterBlock

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

#to embed into comprehensive combustion reactor unit model
# total_flow=1
# biomass_frac=0.1
M_bm = 100
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

m.fs.react.config.reaction_package.rate_reaction_stoichiometry["Rbiomass", "Sol", "uncombustible"].fix(stoich_ash)

#mole_frac_comp spec
m.fs.react.inlet.mole_frac_comp[0,"O2"].fix(N_air*0.21/N_total)
m.fs.react.inlet.mole_frac_comp[0,"N2"].fix(N_air*0.79/N_total)
m.fs.react.inlet.mole_frac_comp[0,"CH4"].fix(1e-20)
m.fs.react.inlet.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.react.inlet.mole_frac_comp[0,"H2O"].fix(1e-20)
m.fs.react.inlet.mole_frac_comp[0,"biomass"].fix(N_bm/N_total)
m.fs.react.inlet.mole_frac_comp[0,"uncombustible"].fix(1e-20)
m.fs.react.inlet.flow_mol.fix(N_total)
m.fs.react.inlet.temperature.fix(30+273.15)
m.fs.react.inlet.pressure.fix(101325)


m.fs.react.rate_reaction_extent[0,"RCH4"].fix(0)
m.fs.react.rate_reaction_extent[0,"Rbiomass"].fix(N_bm/N_total)

# m.fs.react.heat_duty[0].fix()
m.fs.react.outlet.temperature.fix(300+273.15)

print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0
m.fs.react.initialize(outlvl=idaeslog.INFO)
solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)
m.fs.react.report()
