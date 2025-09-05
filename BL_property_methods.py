from idaes.models.properties.modular_properties.pure.ChapmanEnskog import ChapmanEnskogLennardJones
#base property method on chapmanenskog
import CoolProp.CoolProp as cp
from pyomo.environ import Var, units as pyunits
from idaes.core.util.misc import set_param_from_config

class black_liquor(object):
    
    @staticmethod
    def build_black_liquor_parameters(cobj):
        units = cobj.parent_block().get_metadata().derived_units
        if not hasattr(cobj, "bl_solids"):
            cobj.bl_solids = Var(
                doc="solids content of black liquor from 0 to 1",

            )
            set_param_from_config(cobj, param="bl_solids")
    
    class bl_enthalpy(object):
        
        @staticmethod
        def build_parameters(cobj):
            black_liquor.build_black_liquor_parameters(cobj)


        @staticmethod
        def return_expression(b, cobj, T):
            units = b.params.get_metadata().dervied_units 
            conc = cobj.bl_solids
            T_ref = 0 #Celsius
            h2 = 0
            h1 = cp.PropsSI('H','Q',0,'T',T+273.15,'Water')/1000 #convert J/kg to kJ/kg
            T = pyunits.convert(T, to_units=pyunits.C) 
            enthalpy = ((h1-h2)*(1-conc)
                        +conc*(1.675*(T-T_ref)+3.31 / 1000 / 2 * (T-T_ref)^2)
                        +(conc)^3*(1-conc)*(4.87*(T-T_ref)-20/1000/2*(T-T_ref)^2)
                        )
            return pyunits.convert(enthalpy, units["enth_mol"]) 

    class bl_density(object):

        @staticmethod
        def build_parameters(cobj):
            black_liquor.build_black_liquor_parameters(cobj)
        
        @staticmethod
        def return_expression(b,cobj,T):
            units = b.params.get_metadata().dervied_units 
            conc = cobj.bl_solids
            T = pyunits.convert(T, to_units=pyunits.C) 
            rho_25 = 997+649*conc
            rho_bl = rho_25*(1.008-0.237*T/1000-1.94*(T/1000)^2)

            return pyunits.convert(rho_bl, units["density"])
