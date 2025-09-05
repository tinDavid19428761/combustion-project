from pyomo.environ import units as pyunits, Var
# enth_unit = units["flow_mol"]
# flow_unit = Var(1, units=pyunits.m**3/pyunits.s)
# units=[]
flow_unit = units["flow_mol"]
print(f"Flow unit: {pyunits.get_units(flow_unit)}")
