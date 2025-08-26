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
Flash
)
from idaes.models_extra.power_generation.unit_models.helm.phase_separator import HelmPhaseSeparator

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

from custom_combustion_reactor import MultiCombReactor
import unittest

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.flue_properties = GenericParameterBlock(**configuration)

m.fs.steam_properties = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE,
        phase_presentation=PhaseType.LG,
    )

m.fs.R101 = MultiCombReactor(
    property_package = m.fs.flue_properties,
    # reaction_package = m.fs.reaction_params,
    has_heat_of_reaction=True,
    has_heat_transfer=True, 
    has_pressure_change=False,
    has_uncombustibles=True
)

m.fs.R101.conversion["Rbiomass"].fix(0.5)
m.fs.R101.conversion["RCH4"].fix(1)

m.fs.R101.hcon.fix(0.06) #h and w in dh_rxn calculation
# m.fs.R101.wcon.fix(0.09)
# m.fs.R101.wcon.fix(0.11)
m.fs.R101.reaction_package.dh_rxn["Rbiomass"].fix(-2.7804e+06)
# m.fs.R101.ncv.fix(-2749556.40)

m.fs.R101.ash_mass["Rbiomass"].fix(0.03)
# m.fs.R101.ash_mass["RCH4"].fix(0.0)


m.fs.R101.heat_duty[0].fix(-000)
m.fs.R101.surface_area.fix(0.1)
m.fs.R101.surface_temp.fix(55+273.15)

#reactor feed stream
m.fs.R101.inlet.mole_frac_comp[0,"N2"].fix(0.18)
m.fs.R101.inlet.mole_frac_comp[0,"O2"].fix(0.8)
m.fs.R101.inlet.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.R101.inlet.mole_frac_comp[0,"H2O"].fix(1e-20) 
m.fs.R101.inlet.mole_frac_comp[0,"CO"].fix(1e-20) 
m.fs.R101.inlet.mole_frac_comp[0,"biomass"].fix(0.01) 
m.fs.R101.inlet.mole_frac_comp[0,"uncombustible"].fix(1e-20)
m.fs.R101.inlet.mole_frac_comp[0,"CH4"].fix(0.01)
m.fs.R101.inlet.temperature.fix(300)
m.fs.R101.inlet.pressure.fix(101325)
m.fs.R101.inlet.flow_mol.fix(40)

print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0
m.fs.R101.initialize(outlvl=idaeslog.INFO)
solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)
m.fs.R101.report()
