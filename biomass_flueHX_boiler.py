#Importing required pyomo and idaes components
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

#Todo add the four other unit operations
from idaes.models.unit_models import (
Mixer,
StoichiometricReactor,
Heater,
HeatExchanger,
Separator
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
from idaes.models.unit_models.heat_exchanger import delta_temperature_amtd_callback, HX0DInitializer
from idaes.models.unit_models.separator import SplittingType


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

#combustion reactor
m.fs.fire_side = StoichiometricReactor(
    property_package = m.fs.biomass_properties,
    reaction_package = m.fs.reaction_params,
    has_heat_of_reaction=True,
    has_heat_transfer=True,
    has_pressure_change=False,
)

#flue-water heat exchanger
m.fs.boiler_hx = HeatExchanger(
    delta_temperature_callback=delta_temperature_amtd_callback,
    hot_side_name="shell",
    cold_side_name="tube",
    shell={"property_package": m.fs.biomass_properties},
    tube={"property_package": m.fs.steam_properties}
)

m.fs.economizer = HeatExchanger(
    delta_temperature_callback=delta_temperature_amtd_callback,
    hot_side_name="shell",
    cold_side_name="tube",
    shell={"property_package": m.fs.biomass_properties},
    tube={"property_package": m.fs.steam_properties}
)

m.fs.ash_sep = Separator(
    property_package = m.fs.biomass_properties,
    split_basis = SplittingType.phaseFlow,
    outlet_list = ["flue", "ash"]
)

m.fs.bdw_sep = Separator(
    property_package = m.fs.steam_properties,
    split_basis = SplittingType.totalFlow,
    ideal_separation = False,
    outlet_list = ["bdw", "to_boiler"]
)

#reactor flow sheet: feed-> reactor -> hot flue -> HX -> flue
m.fs.s01 = Arc(source=m.fs.fire_side.outlet,destination=m.fs.ash_sep.inlet)
m.fs.s02 = Arc(source=m.fs.ash_sep.flue,destination=m.fs.boiler_hx.shell_inlet)
m.fs.s03 = Arc(source=m.fs.boiler_hx.shell_outlet,destination=m.fs.economizer.shell_inlet)
m.fs.s04 = Arc(source=m.fs.economizer.tube_outlet,destination=m.fs.bdw_sep.inlet)
m.fs.s05 = Arc(source=m.fs.bdw_sep.to_boiler,destination=m.fs.boiler_hx.tube_inlet)
TransformationFactory("network.expand_arcs").apply_to(m)

#specifying reaction package parameters
m.fs.reaction_params.h.fix(0.06) #h and w in dh_rxn calculation
m.fs.reaction_params.w.fix(0.09)
m.fs.reaction_params.rate_reaction_stoichiometry["R1","Sol","ash"]=0.02 #specify ash content

#reaction conversion constraint
m.fs.fire_side.conversion = Var(initialize=1, bounds=(0,1))
m.fs.fire_side.conversion_constraint = Constraint(
    expr=m.fs.fire_side.conversion*m.fs.fire_side.inlet.mole_frac_comp[0,"biomass"]*m.fs.fire_side.inlet.flow_mol[0]
    == (
        m.fs.fire_side.inlet.flow_mol[0]*m.fs.fire_side.inlet.mole_frac_comp[0,"biomass"]
        -m.fs.fire_side.outlet.flow_mol[0]*m.fs.fire_side.outlet.mole_frac_comp[0,"biomass"]
    )
)
m.fs.fire_side.conversion.fix(1)

#modelling Q_loss_casing_convection
m.fs.fire_side.ohtc = Param(initialize=250, units=pyunits.J/pyunits.m**2/pyunits.K/pyunits.s, doc="boiler casing overall heat transfer coefficient")
m.fs.fire_side.surface_area = Param(initialize=0.02, units=pyunits.m**2, doc="casing outer surface area")
m.fs.fire_side.surface_temp = Param(initialize=55+273.12, units=pyunits.K, doc="outer skin temperature of boiler")

m.fs.fire_side.heat_loss_constraint = Constraint(
    expr=-m.fs.fire_side.heat_duty[0]==
    (m.fs.fire_side.ohtc*m.fs.fire_side.surface_area*(m.fs.fire_side.outlet.temperature[0]-m.fs.fire_side.surface_temp))
)

#reactor feed stream
m.fs.fire_side.inlet.mole_frac_comp[0,"N2"].fix(0.7)
m.fs.fire_side.inlet.mole_frac_comp[0,"O2"].fix(0.29)
m.fs.fire_side.inlet.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.fire_side.inlet.mole_frac_comp[0,"H2O"].fix(1e-20) 
m.fs.fire_side.inlet.mole_frac_comp[0,"CO"].fix(1e-20) 
m.fs.fire_side.inlet.mole_frac_comp[0,"biomass"].fix(0.01) 
m.fs.fire_side.inlet.mole_frac_comp[0,"ash"].fix(1e-20)
m.fs.fire_side.inlet.temperature.fix(400)
m.fs.fire_side.inlet.pressure.fix(101325)
m.fs.fire_side.inlet.flow_mol.fix(40)

#specifying ash separation
m.fs.ash_sep.split_fraction[0,"ash","Sol"].fix(1)
m.fs.ash_sep.split_fraction[0,"ash","Vap"].fix(0)

m.fs.economizer.tube_inlet.enth_mol.fix(m.fs.steam_properties.htpx(p=101325*pyunits.Pa,T=350*pyunits.K))
m.fs.economizer.tube_outlet.enth_mol.fix(m.fs.steam_properties.htpx(p=101325*pyunits.Pa,x=0))
m.fs.boiler_hx.tube_outlet.enth_mol.fix(m.fs.steam_properties.htpx(p=101325*pyunits.Pa,T=400*pyunits.K))
m.fs.boiler_hx.overall_heat_transfer_coefficient[0].fix(100)
m.fs.economizer.tube_inlet.flow_mol.fix(20)
m.fs.economizer.tube_inlet.pressure.fix(101325)
m.fs.economizer.overall_heat_transfer_coefficient[0].fix(100)

m.fs.bdw_sep.split_fraction[0,"bdw"].fix(0.05)

#initialization routine
seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"
seq.options.tear_method = "Wegstein"
seq.options.iterLim = 3

G = seq.create_graph(m)
heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
order = seq.calculation_order(G)

#for identifying tear stream
""" for o in heuristic_tear_set: #fs.s03
    print(o.name) """ 

tear_guesses = {
    "mole_frac_comp": {
        (0, "N2"): 0.67,
        (0, "O2"): 0.22,
        (0, "CO2"): 0.06,
        (0, "H2O"): 0.05,
        (0, "CO"): 1e-20,
        (0, "biomass"): 1e-20,
        (0, "ash"): 1e-20,
    },
    "flow_mol": {0: 40},
    "temperature": {0: 610},
    "pressure": {0: 101325},
}

seq.set_guesses_for(m.fs.economizer.shell_inlet, tear_guesses)

def function(unit):
    unit.initialize(outlvl=idaeslog.INFO)

print(degrees_of_freedom(m))  
seq.run(m, function)

#pre-solve [actual] re-specification
m.fs.fire_side.inlet.flow_mol.unfix()
m.fs.economizer.shell_outlet.temperature.fix(500)

solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)

print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0

m.fs.boiler_eff = Expression(
    expr = (m.fs.boiler_hx.heat_duty[0]+m.fs.economizer.heat_duty[0])/(m.fs.fire_side.rate_reaction_extent[0,"R1"]*-m.fs.reaction_params.dh_rxn["R1"])
)

#results
m.fs.fire_side.report()
m.fs.boiler_hx.report()
m.fs.economizer.report()
print(f"    Boiler Efficiency: {value(m.fs.boiler_eff)*100:.2f}%")
print(f"    direct heat loss:{value(m.fs.fire_side.heat_duty[0]):.2f}J/s")



