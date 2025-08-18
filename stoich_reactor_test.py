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

from custom_stoichiometric_reactor import StoichiometricReactor
import unittest

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.biomass_properties = GenericParameterBlock(**configuration)
# m.fs.reaction_params = BMCombReactionParameterBlock(property_package=m.fs.biomass_properties)

m.fs.steam_properties = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE,
        phase_presentation=PhaseType.LG,
        # state_vars=StateVars.TPX
    )

m.fs.R101 = StoichiometricReactor(
    property_package = m.fs.biomass_properties,
    # reaction_package = m.fs.reaction_params,
    has_heat_of_reaction=True,
    has_heat_transfer=False, #test with true also
    has_pressure_change=False,
)

m.fs.R101.conversion.fix(0.5)
m.fs.R101.reaction_package.h.fix(0.06) #h and w in dh_rxn calculation
m.fs.R101.reaction_package.w.fix(0.09)

m.fs.R101.reaction_package.rate_reaction_stoichiometry["R1","Sol","ash"].fix(0.03)


#reactor feed stream
m.fs.R101.inlet.mole_frac_comp[0,"N2"].fix(0.7)
m.fs.R101.inlet.mole_frac_comp[0,"O2"].fix(0.29)
m.fs.R101.inlet.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.R101.inlet.mole_frac_comp[0,"H2O"].fix(1e-20) 
m.fs.R101.inlet.mole_frac_comp[0,"CO"].fix(1e-20) 
m.fs.R101.inlet.mole_frac_comp[0,"biomass"].fix(0.01) 
m.fs.R101.inlet.mole_frac_comp[0,"ash"].fix(1e-20)
m.fs.R101.inlet.temperature.fix(300)
m.fs.R101.inlet.pressure.fix(101325)
m.fs.R101.inlet.flow_mol.fix(40)

print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0
m.fs.R101.initialize(outlvl=idaeslog.INFO)
solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)
m.fs.R101.report()
print(value(m.fs.R101.reaction_package.rate_reaction_stoichiometry["R1", "Sol", "ash"]))
print(value(m.fs.R101.reaction_package.dh_rxn["R1"]))

# assert value(m.fs.R101.reaction_package.dh_rxn["R1"]) == approx(-2749556.4, rel=1e-6)
# assert value(m.fs.R101.outlet.temperature) == approx(727.15, rel=1e-3)