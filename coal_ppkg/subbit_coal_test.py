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
StoichiometricReactor
)
from idaes.models_extra.power_generation.unit_models.helm.phase_separator import HelmPhaseSeparator

from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
# Import idaes logger to set output levels
import idaes.logger as idaeslog
from idaes.models.properties.modular_properties import GenericParameterBlock, GenericReactionParameterBlock


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
from idaes.core.util.model_diagnostics import (DiagnosticsToolbox,)

import unittest

from subbit_coal_ppkg import configuration 
from subbit_coal_rxnpkg import rxn_configuration

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.coal_properties = GenericParameterBlock(**configuration)
m.fs.reaction_params = GenericReactionParameterBlock(property_package=m.fs.coal_properties, **rxn_configuration)

m.fs.R101 = StoichiometricReactor(
    property_package = m.fs.coal_properties,
    reaction_package = m.fs.reaction_params,
    has_heat_of_reaction=True,
    has_heat_transfer=True,
    has_pressure_change=False,
)

m.fs.R101.rate_reaction_extent[0,"Rcoal"].fix(0.1) #half of 0.1 mol/s of coal


#reactor feed stream
m.fs.R101.inlet.mole_frac_comp[0,"coal"].fix(0.05) 
m.fs.R101.inlet.mole_frac_comp[0,"N2"].fix(0.65)
m.fs.R101.inlet.mole_frac_comp[0,"O2"].fix(0.2)
m.fs.R101.inlet.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.R101.inlet.mole_frac_comp[0,"NO"].fix(1e-20)
m.fs.R101.inlet.mole_frac_comp[0,"H2O"].fix(1e-20) 
m.fs.R101.inlet.mole_frac_comp[0,"uncombustible"].fix(1e-20)
m.fs.R101.inlet.temperature.fix(300)
m.fs.R101.inlet.pressure.fix(101325)
m.fs.R101.inlet.flow_mol.fix(10)

m.fs.R101.outlet.temperature.fix(300)
# m.fs.R101.heat_duty.fix(0)

print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0
m.fs.R101.initialize(outlvl=idaeslog.INFO)
solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)
m.fs.R101.report()




