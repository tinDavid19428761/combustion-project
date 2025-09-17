# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 14:18:30 2023

@author: bjl25

validation notes:
validated in excel
heat of reaction validated as heat of formations
temperature raise and boiling validated
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

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.biomass_properties = GenericParameterBlock(**configuration)
m.fs.reaction_params = BMCombReactionParameterBlock(
    property_package=m.fs.biomass_properties
)

m.fs.steam_properties = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE,
        phase_presentation=PhaseType.LG,
    )

from idaes.core.util.model_statistics import degrees_of_freedom


# creating reactor unit
m.fs.M101 = Mixer(
    property_package = m.fs.biomass_properties,
    inlet_list=["biomass_feed"]
)

m.fs.R101 = StoichiometricReactor(
    property_package = m.fs.biomass_properties,
    reaction_package = m.fs.reaction_params,
    has_heat_of_reaction=True,
    has_heat_transfer=False,
    has_pressure_change=False,
)

m.fs.H101 = Heater(
    property_package = m.fs.steam_properties,
    has_pressure_change=False,
    has_phase_equilibrium=False,
)


#reactor flow sheet feed via mixer -> reactor -> product via separator

m.fs.s01 = Arc(source=m.fs.M101.outlet,destination=m.fs.R101.inlet)

TransformationFactory("network.expand_arcs").apply_to(m)


flowTotal = 1
extentR1 = 0.01*flowTotal #absolute extent of reaction

m.fs.M101.biomass_feed.mole_frac_comp[0,"N2"].fix(0.5)
m.fs.M101.biomass_feed.mole_frac_comp[0,"O2"].fix(0.47)
m.fs.M101.biomass_feed.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.M101.biomass_feed.mole_frac_comp[0,"H2O"].fix(1e-20) 
m.fs.M101.biomass_feed.mole_frac_comp[0,"CO"].fix(1e-20) 
m.fs.M101.biomass_feed.mole_frac_comp[0,"biomass"].fix(0.03) 
m.fs.M101.biomass_feed.temperature.fix(300)
m.fs.M101.biomass_feed.pressure.fix(101325)
m.fs.M101.biomass_feed.flow_mol.fix(flowTotal)

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
m.fs.R101.outlet.temperature.unfix()

m.fs.R101.initialize()
#solve reactor
solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)

#steam generation specification
steam_pressure = 101325
T_water_in=300
T_water_out=380
enth_in=m.fs.steam_properties.htpx(p=steam_pressure * pyunits.Pa, T= T_water_in * pyunits.K)
enth_out=m.fs.steam_properties.htpx(p=steam_pressure * pyunits.Pa, T= T_water_out * pyunits.K)

#steam_gen for initialization purposes
initSteam = 10

m.fs.H101.inlet.flow_mol.fix(initSteam)
m.fs.H101.inlet.enth_mol.fix(enth_in)
m.fs.H101.inlet.pressure.fix(steam_pressure)
# m.fs.H101.heat_duty.fix(-value(m.fs.R101.heat_duty[0]))
m.fs.H101.heat_duty.fix(12345)

m.fs.H101.initialize(outlvl=idaeslog.INFO)

#re-specifying .fix's
m.fs.H101.inlet.flow_mol.unfix()
m.fs.H101.outlet.enth_mol.fix(enth_out)

#solve fow steam generation flow
solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)
# print(value(m.fs.R101.heat_duty[0]))
print(value(m.fs.H101.heat_duty[0]))

m.fs.R101.report()
m.fs.H101.report()
