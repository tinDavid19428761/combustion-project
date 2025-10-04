import CoolProp.CoolProp as CP
fluid='methane'
enth = CP.PropsSI('H','T',273.15+50,'P',101325,fluid)
print(enth)