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

from biomass_prefix.biomass_combustion_rp import BMCombReactionParameterBlock

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

from biomass_prefix.biomass_comb_pp import configuration
from modular_biomass_rxn import rxn_configuration

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.properties = GenericParameterBlock(**configuration)
m.fs.reaction = GenericReactionParameterBlock(property_package=m.fs.properties, **rxn_configuration)


m.fs.react = StoichiometricReactor(
    property_package = m.fs.properties,
    reaction_package = m.fs.reaction,
    has_heat_of_reaction=True,
    has_heat_transfer=True, 
    has_pressure_change=False
)

m.fs.react.inlet.mole_frac_comp[0,"O2"].fix(0.28)
m.fs.react.inlet.mole_frac_comp[0,"N2"].fix(0.7)
m.fs.react.inlet.mole_frac_comp[0,"biomass"].fix(0.02)
m.fs.react.inlet.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.react.inlet.mole_frac_comp[0,"H2O"].fix(1e-20)
m.fs.react.inlet.mole_frac_comp[0,"uncombustible"].fix(1e-20)
m.fs.react.inlet.flow_mol.fix(10)
m.fs.react.inlet.temperature.fix(30+273.15)
m.fs.react.inlet.pressure.fix(101325)

m.fs.react.rate_reaction_extent[0,"Rbiomass"].fix(0.2)

# m.fs.react.heat_duty[0].fix()
m.fs.react.outlet.temperature.fix(300+273.15)

print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0
m.fs.react.initialize(outlvl=idaeslog.INFO)
solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)
m.fs.react.report()
