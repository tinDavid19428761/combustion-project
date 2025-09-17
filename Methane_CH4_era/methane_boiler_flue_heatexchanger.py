# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 14:18:30 2023

@author: bjl25
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
HeatExchanger
)

from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
# Import idaes logger to set output levels
import idaes.logger as idaeslog
from idaes.models.properties.modular_properties import GenericParameterBlock
from  CH4_Comb_PPwithCO import configuration 
from  ch4_comb_reactionWithoutKinetics import CH4CombReactionParameterBlock

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

m.fs.methane_properties = GenericParameterBlock(**configuration)
m.fs.reaction_params = CH4CombReactionParameterBlock(
    property_package=m.fs.methane_properties
)

m.fs.steam_properties = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE,
        phase_presentation=PhaseType.LG,
    )



# creating reactor unit
m.fs.M101 = Mixer(
    property_package = m.fs.methane_properties,
    inlet_list=["methane_feed"]
)

m.fs.R101 = StoichiometricReactor(
    property_package = m.fs.methane_properties,
    reaction_package = m.fs.reaction_params,
    has_heat_of_reaction=True,
    has_heat_transfer=True,
    has_pressure_change=False,
)

#flue-water heat exchanger
m.fs.E101 = HeatExchanger(
    delta_temperature_callback=delta_temperature_amtd_callback,
    hot_side_name="shell",
    cold_side_name="tube",
    shell={"property_package": m.fs.methane_properties},
    tube={"property_package": m.fs.steam_properties}
)

#reactor flow sheet feed via mixer -> reactor -> product via separator

m.fs.s01 = Arc(source=m.fs.M101.outlet,destination=m.fs.R101.inlet)
m.fs.s02 = Arc(source=m.fs.R101.outlet,destination=m.fs.E101.shell_inlet)

TransformationFactory("network.expand_arcs").apply_to(m)


def O2feed(flowCH4,excess):
    return flowCH4*2*(1+excess), flowCH4*2*79/21*(1+excess)

# check sufficient O2 for given feed and completeness, then calc R1,R2 reaction extents
def stoichExtents(flowCH4,flowO2,complete): #returns extentCO2, extentCO

    if flowCH4*((complete*2)+((1-complete)*1.5))<=flowO2:
        return complete*flowCH4, (1-complete)*flowCH4
    else:
        return complete*(flowO2/((complete*2)+((1-complete)*1.5))), (1-complete)*(flowO2/((complete*2)+((1-complete)*1.5)))


flowCH4=0.2
completeness=0.7
excessAir=0.5

flowO2,flowN2 = O2feed(flowCH4,excessAir)
extentCO2, extentCO = stoichExtents(flowCH4,flowO2,completeness)

flowTotal=sum([flowO2,flowCH4,flowN2])

m.fs.M101.methane_feed.mole_frac_comp[0,"N2"].fix(flowN2/flowTotal)
m.fs.M101.methane_feed.mole_frac_comp[0,"O2"].fix(flowO2/flowTotal)
m.fs.M101.methane_feed.mole_frac_comp[0,"CH4"].fix(flowCH4/flowTotal)
m.fs.M101.methane_feed.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.M101.methane_feed.mole_frac_comp[0,"H2O"].fix(1e-20) 
m.fs.M101.methane_feed.mole_frac_comp[0,"CO"].fix(1e-20) 
m.fs.M101.methane_feed.temperature.fix(1000)
m.fs.M101.methane_feed.pressure.fix(10000)
m.fs.M101.methane_feed.flow_mol.fix(flowTotal)


# following ~20 lines of code courtesy of chatgpt
# Add variables for extent of each reaction (mol/s)
m.fs.R101.extent_R1 = Var(m.fs.time, initialize=extentCO2, units=pyunits.mol/pyunits.s)
m.fs.R101.extent_R2 = Var(m.fs.time, initialize=extentCO, units=pyunits.mol/pyunits.s)

# Override the default extent expressions
def extent_match_rule(b, t, r):
    if r == "R1":
        return b.rate_reaction_extent[t, r] == b.extent_R1[t]
    elif r == "R2":
        return b.rate_reaction_extent[t, r] == b.extent_R2[t]
    else:
        return Constraint.Skip

m.fs.R101.extent_match = Constraint(
    m.fs.time, m.fs.R101.config.reaction_package.rate_reaction_idx,
    rule=extent_match_rule
)

m.fs.R101.extent_R1[0.0].fix(extentCO2) 
m.fs.R101.extent_R2[0.0].fix(extentCO) 
m.fs.R101.outlet.temperature.unfix()
m.fs.R101.heat_duty.fix(0)

print(degrees_of_freedom(m))

m.fs.R101.initialize()

#specifying heat exchanger
m.fs.E101.area.fix(0.1)
m.fs.E101.overall_heat_transfer_coefficient[0].fix(100)
m.fs.E101.tube_inlet.flow_mol.fix(0.2)
m.fs.E101.tube_inlet.pressure.fix(101325)
m.fs.E101.tube_inlet.enth_mol.fix(m.fs.steam_properties.htpx(p=101325*pyunits.Pa,T=290*pyunits.K))

initializer = HX0DInitializer()
initializer.initialize(m.fs.E101)

# re-specifying. Solves for 
m.fs.E101.area.unfix()
m.fs.E101.tube_outlet.enth_mol.fix(m.fs.steam_properties.htpx(p=101325*pyunits.Pa,T=400*pyunits.K))
m.fs.E101.tube_inlet.flow_mol.unfix()
m.fs.E101.shell_outlet.temperature.fix(1100)
print(degrees_of_freedom(m))

solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)

print(degrees_of_freedom(m))

# print(value(m.fs.R101.outlet.temperature[0]))
m.fs.E101.report()
m.fs.R101.report()