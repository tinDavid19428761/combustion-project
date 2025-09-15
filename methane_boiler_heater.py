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
from  CH4_Comb_PPwithCO import configuration 
from  ch4_comb_reactionWithoutKinetics import CH4CombReactionParameterBlock

#helmholtz import for water
from idaes.models.properties.general_helmholtz import (
        HelmholtzParameterBlock,
        HelmholtzThermoExpressions,
        AmountBasis,
        PhaseType,
    )

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

# m.fs.dsi.properties_steam_in[0].enth_mol.fix(
#     m.fs.steam_properties.htpx(p=101325 * pyunits.Pa, T= 300 * pyunits.K)
# )

from idaes.core.util.model_statistics import degrees_of_freedom


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

m.fs.H101 = Heater(
    property_package = m.fs.steam_properties,
    has_pressure_change=False,
    has_phase_equilibrium=False,
)


#reactor flow sheet feed via mixer -> reactor -> product via separator

m.fs.s01 = Arc(source=m.fs.M101.outlet,destination=m.fs.R101.inlet)

TransformationFactory("network.expand_arcs").apply_to(m)

# check sufficient O2 for given feed and completeness, then calc R1,R2 reaction extents
def stoichExtents(flowCH4,flowO2,complete): #returns extentCO2, extentCO

    if flowCH4*((complete*2)+((1-complete)*1.5))<=flowO2:
        return complete*flowCH4, (1-complete)*flowCH4
    else:
        return complete*(flowO2/((complete*2)+((1-complete)*1.5))), (1-complete)*(flowO2/((complete*2)+((1-complete)*1.5)))

def O2feed(flowCH4,excess):
    return flowCH4*2*(1+excess), flowCH4*2*79/21*(1+excess)

flowCH4=0.2
completeness=0.7
excessAir=0.5

#hard-specified air feed.
# flowO2=0.3
# flowN2=flowO2*79/21

#air feed specified by excess air:
flowO2,flowN2 = O2feed(flowCH4,excessAir)
extentCO2, extentCO = stoichExtents(flowCH4,flowO2,completeness)

flowTotal=sum([flowO2,flowCH4,flowN2])

m.fs.M101.methane_feed.mole_frac_comp[0,"N2"].fix(flowN2/flowTotal)
m.fs.M101.methane_feed.mole_frac_comp[0,"O2"].fix(flowO2/flowTotal)
m.fs.M101.methane_feed.mole_frac_comp[0,"CH4"].fix(flowCH4/flowTotal)
m.fs.M101.methane_feed.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.M101.methane_feed.mole_frac_comp[0,"H2O"].fix(1e-20) 
m.fs.M101.methane_feed.mole_frac_comp[0,"CO"].fix(1e-20) 
m.fs.M101.methane_feed.temperature.fix(1100)
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
m.fs.R101.outlet.temperature.fix(1099.9)

m.fs.R101.initialize(outlvl=idaeslog.INFO)
#solve reactor
solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)

#_______________________

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
m.fs.H101.heat_duty.fix(-value(m.fs.R101.heat_duty[0]))

m.fs.H101.initialize(outlvl=idaeslog.INFO)

#re-specifying .fix's
m.fs.H101.inlet.flow_mol.unfix()
m.fs.H101.outlet.enth_mol.fix(enth_out)
#solve fow steam generation flow
solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)

m.fs.R101.report()
m.fs.H101.report()
