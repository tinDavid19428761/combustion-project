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
HeatExchanger
)

from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
# Import idaes logger to set output levels
import idaes.logger as idaeslog
from idaes.models.properties.modular_properties import GenericParameterBlock
from  single_comp_biomass_comb_pp import configuration 
from  biomass_combustion_rp import BMCombReactionParameterBlock

#helmholtz import for water
from idaes.models.properties.general_helmholtz import (
        HelmholtzParameterBlock,
        HelmholtzThermoExpressions,
        AmountBasis,
        PhaseType,
    )
from idaes.models.unit_models.heat_exchanger import delta_temperature_amtd_callback, HX0DInitializer

from idaes.core.util.model_statistics import degrees_of_freedom

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.biomass_properties = GenericParameterBlock(**configuration)
m.fs.reaction_params = BMCombReactionParameterBlock(
    property_package=m.fs.biomass_properties
)
m.fs.reaction_params.h.fix(0.06)
m.fs.reaction_params.w.fix(0.09)

m.fs.steam_properties = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE,
        phase_presentation=PhaseType.LG,
    )

#combustion reactor
m.fs.R101 = StoichiometricReactor(
    property_package = m.fs.biomass_properties,
    reaction_package = m.fs.reaction_params,
    has_heat_of_reaction=True,
    has_heat_transfer=False,
    has_pressure_change=False,
)

#flue-water heat exchanger
m.fs.E101 = HeatExchanger(
    delta_temperature_callback=delta_temperature_amtd_callback,
    hot_side_name="shell",
    cold_side_name="tube",
    shell={"property_package": m.fs.biomass_properties},
    tube={"property_package": m.fs.steam_properties}
)

#reactor flow sheet: feed-> reactor -> hot flue -> HX -> flue

m.fs.s02 = Arc(source=m.fs.R101.outlet,destination=m.fs.E101.shell_inlet)
TransformationFactory("network.expand_arcs").apply_to(m)

#specifying reactor
flowTotal = 1
extentR1 = 0.01*flowTotal #absolute extent of reaction

m.fs.R101.inlet.mole_frac_comp[0,"N2"].fix(0.7)
m.fs.R101.inlet.mole_frac_comp[0,"O2"].fix(0.29)
m.fs.R101.inlet.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.R101.inlet.mole_frac_comp[0,"H2O"].fix(1e-20) 
m.fs.R101.inlet.mole_frac_comp[0,"CO"].fix(1e-20) 
m.fs.R101.inlet.mole_frac_comp[0,"biomass"].fix(0.01) 
m.fs.R101.inlet.temperature.fix(500)
m.fs.R101.inlet.pressure.fix(101325)
m.fs.R101.inlet.flow_mol.fix(flowTotal)

# Add variables for extent of each reaction (mol/s)
m.fs.R101.extent_R1 = Var(m.fs.time, initialize=extentR1, units=pyunits.mol/pyunits.s)

# Override the default extent expressions
def extent_match_rule(b, t, r):
    if r == "R1":
        return b.rate_reaction_extent[t, r] == b.extent_R1[t]
    else:
        return Constraint.Skip

m.fs.R101.extent_match = Constraint(
    m.fs.time, m.fs.R101.config.reaction_package.rate_reaction_idx,
    rule=extent_match_rule
)

m.fs.R101.extent_R1[0.0].fix(extentR1)

#initialisation routine:
#heat exchanger init specs
m.fs.E101.area.fix(0.5)
m.fs.E101.overall_heat_transfer_coefficient[0].fix(100)
m.fs.E101.tube_inlet.flow_mol.fix(0.5)
m.fs.E101.tube_inlet.pressure.fix(101325)
m.fs.E101.tube_inlet.enth_mol.fix(m.fs.steam_properties.htpx(p=101325*pyunits.Pa,T=290*pyunits.K))
#actual init.
seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"

def function(unit):
    unit.initialize(outlvl=idaeslog.INFO)
    
seq.run(m, function)

#final heat exchanger specifying
m.fs.E101.shell_outlet.temperature.fix(600)
# m.fs.E101.area.unfix()
m.fs.E101.overall_heat_transfer_coefficient[0].unfix()
m.fs.E101.tube_inlet.flow_mol.unfix()
m.fs.E101.tube_outlet.enth_mol.fix(m.fs.steam_properties.htpx(p=101325*pyunits.Pa,T=400*pyunits.K))

print(degrees_of_freedom(m))

solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)

print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0

#results
m.fs.E101.report()
m.fs.R101.report()
print(value(m.fs.reaction_params.ncv))
print(value(m.fs.R101.extent_R1[0.0]))