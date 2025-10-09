#black liquor case study

import matplotlib.pyplot as plt

#stoich_reactor_test
from pyomo.environ import (
    Constraint,
    Var,
    ConcreteModel,
    Expression,
    SolverFactory,
    TransformationFactory,
    value,
    units as pyunits
)

#Todo add the four other unit operations
from idaes.models.unit_models import (
Mixer,
# StoichiometricReactor,
Heater,
Heater
)

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog
from idaes.models.properties.modular_properties import GenericParameterBlock
from black_liquor_pp import configuration

#helmholtz import for water
from idaes.models.properties.general_helmholtz import (
        HelmholtzParameterBlock,
        AmountBasis,
        PhaseType,
    )

from bl_combustion_reactor import MultiCombReactor #has rxn pkg included

from pyomo.environ import Reference, Var, Param, units as pyunits, value

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.bl_properties = GenericParameterBlock(**configuration)

m.fs.steam_properties = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE,
        phase_presentation=PhaseType.LG,
    )

m.fs.R101 = MultiCombReactor(
    property_package = m.fs.bl_properties,
    # reaction_package = m.fs.reaction_params,
    has_heat_of_reaction=True,
    has_heat_transfer=True, 
    has_pressure_change=False
)

m.fs.H101 = Heater(
    property_package = m.fs.steam_properties,
)

m.fs.mix = Mixer(
    property_package = m.fs.bl_properties,
    inlet_list = ["air","fuel"]
)

m.fs.s01 = Arc(source=m.fs.mix.outlet, destination=m.fs.R101.inlet)
TransformationFactory("network.expand_arcs").apply_to(m)


#case study inputs
#flows
f = {
    "air_flow": 976,
    "bl_flow": 998,
    "gas_flow": 0.1717
}


fuel_total = f["bl_flow"]+f["gas_flow"]


# air stream
m.fs.mix.air.mole_frac_comp[0,"N2"].fix(0.79)
m.fs.mix.air.mole_frac_comp[0,"O2"].fix(0.21)
m.fs.mix.air.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.mix.air.mole_frac_comp[0,"H2O"].fix(1e-20)
m.fs.mix.air.mole_frac_comp[0,"BL"].fix(1e-20)
m.fs.mix.air.mole_frac_comp[0,"uncombustible"].fix(1e-20)
m.fs.mix.air.mole_frac_comp[0,"CH4"].fix(1e-20)
m.fs.mix.air.temperature.fix(33+273.15)
m.fs.mix.air.pressure.fix(101325)
m.fs.mix.air.flow_mol.fix(f["air_flow"]) 

# black liquor stream stream
m.fs.mix.fuel.mole_frac_comp[0,"N2"].fix(1e-20)
m.fs.mix.fuel.mole_frac_comp[0,"O2"].fix(1e-20)
m.fs.mix.fuel.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.mix.fuel.mole_frac_comp[0,"H2O"].fix(1e-20)
m.fs.mix.fuel.mole_frac_comp[0,"BL"].fix(f["bl_flow"]/fuel_total)
m.fs.mix.fuel.mole_frac_comp[0,"uncombustible"].fix(1e-20)
m.fs.mix.fuel.mole_frac_comp[0,"CH4"].fix(f["gas_flow"]/fuel_total)
m.fs.mix.fuel.temperature.fix(123.6+273.15)
m.fs.mix.fuel.pressure.fix(101325)
m.fs.mix.fuel.flow_mol.fix(fuel_total)

m.fs.R101.conversion_Rbl.fix(1)
m.fs.R101.conversion_RCH4.fix(0)

m.fs.R101.dh_rxn_RCH4.fix(-802125)
m.fs.R101.dh_rxn_Rbl.fix(-135150)

m.fs.R101.outlet.temperature.fix(195.66+273.15)


m.fs.mix.initialize(outlvl=idaeslog.INFO)
m.fs.R101.initialize(outlvl=idaeslog.INFO)

m.fs.H101.inlet.flow_mol.fix(2177.69)
m.fs.H101.inlet.enth_mol.fix(m.fs.steam_properties.htpx(p=45*100*1000*pyunits.Pa,T=(136.39+273.15)*pyunits.K))
m.fs.H101.inlet.pressure.fix(45*100*1000)

m.fs.H101.outlet.enth_mol.fix(m.fs.steam_properties.htpx(p=45*100*1000*pyunits.Pa,T=(400+273.15)*pyunits.K))


m.fs.R101.heat_loss = Var(initialize=100000, units=pyunits.J/pyunits.s)


def heat_transfer_rule(b,t):
    return b.H101.heat_duty[0]+b.R101.heat_loss == -b.R101.heat_duty[0]

m.fs.heat_transfer = Constraint(
    m.fs.time,
    rule=heat_transfer_rule
)


print(degrees_of_freedom(m.fs.H101))
m.fs.H101.initialize(outlvl=idaeslog.INFO)


m.fs.boiler_efficiency = Expression(
        expr = 100*value(m.fs.H101.heat_duty[0])/sum(m.fs.R101.control_volume.rate_reaction_extent[0,r]*(-m.fs.R101.reaction_package.dh_rxn[r]) for r in m.fs.R101.reaction_package.rate_reaction_idx)
    )

solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)
case_fluetemp = value(m.fs.R101.outlet.temperature[0])
case_eff = value(m.fs.boiler_efficiency)

m.fs.R101.report()
print(case_eff)
print(case_fluetemp)
print(value(m.fs.R101.heat_loss))

m.fs.R101.heat_loss.fix(30604244.11776122) #solved for case study condition

m.fs.H101.inlet.flow_mol.unfix()

print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0
temps = range(600-273,150,-10)
effs = [0]*len(temps)


for i,flue_temp in enumerate(temps):
    m.fs.R101.outlet.temperature[0].fix(flue_temp+273.15)
    
    solver=SolverFactory("ipopt")
    status=solver.solve(m,tee=True)

    m.fs.boiler_efficiency = Expression(
        expr = 100*value(m.fs.H101.heat_duty[0])/sum(m.fs.R101.control_volume.rate_reaction_extent[0,r]*(-m.fs.R101.reaction_package.dh_rxn[r]) for r in m.fs.R101.reaction_package.rate_reaction_idx)
    )
    effs[i] = value(m.fs.boiler_efficiency)


m.fs.R101.report()
m.fs.H101.report()
# m.fs.mix.report()
print(f"{value(m.fs.boiler_efficiency):.2f}%")
print(f"{(value(m.fs.R101.outlet.temperature[0])-273.15):.2f} C")
# m.fs.R101.display()

print(case_eff)
print(case_fluetemp)

plt.plot(temps,effs)
plt.ylim(0,100)
plt.xlabel("Flue Stack Temperature [C]")
plt.ylabel("Boiler Efficiency [%]")
plt.title("Boiler Efficiency Curve Against Stack Temperature")
plt.plot([case_fluetemp-273.15],[case_eff],'ro',label='Case Study Conditions')
plt.legend(loc='center left')
plt.show()

