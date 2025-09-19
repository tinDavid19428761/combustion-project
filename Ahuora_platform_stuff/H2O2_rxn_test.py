#H2O2 CSTr test

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
CSTR,
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

from H2O2_ppkg import configuration
from H2O2_first_order_rxnpkg import config_dict

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

m.fs.h2o2_properties = GenericParameterBlock(**configuration)
m.fs.h2o2_rxn_prop = GenericReactionParameterBlock(property_package=m.fs.h2o2_properties, **config_dict)

m.fs.react = CSTR(
    property_package = m.fs.h2o2_properties,
    reaction_package = m.fs.h2o2_rxn_prop,
    has_heat_transfer = True,

)

m.fs.react.volume.fix(1000)

m.fs.react.inlet.mole_frac_comp[0,"H2O2"].fix(0.99)
m.fs.react.inlet.mole_frac_comp[0,"H2O"].fix(0.01)
m.fs.react.inlet.mole_frac_comp[0,"O2"].fix(1e-20)
m.fs.react.inlet.temperature.fix(26+273.15)
m.fs.react.inlet.pressure.fix(101325)
m.fs.react.inlet.flow_mol.fix(10)

m.fs.react.outlet.temperature.fix(50+273.15)

print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0
m.fs.react.initialize(outlvl=idaeslog.INFO)
solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)
m.fs.react.report()
