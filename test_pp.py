"""
For testing viscosity in pp
"""
# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import VaporPhase, Component, SolidPhase
import idaes.logger as idaeslog

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.pure import ConstantProperties, ChapmanEnskogLennardJones 
from idaes.models.properties.modular_properties.pure import ChapmanEnskog
from idaes.models.properties.modular_properties.pure.ChapmanEnskog import collision_integral_neufeld_callback 
from idaes.models.properties.modular_properties.pure import NIST
from idaes.core import PhaseType as PT
from idaes.models_extra.power_generation.properties.flue_gas_ideal import FlueGasStateBlock


# Set up logger
_log = idaeslog.getLogger(__name__)

configuration = {
    # Specifying components
    "components": {
        
        "ash": {
            "type": Component,
            "elemental_composition": {"Random":1}, #mainly SiO2, Al2O3, CaO, Fe2O3, MgO
            "cp_mol_sol_comp": ConstantProperties.Constant,
            "enth_mol_sol_comp": ConstantProperties.Constant,
            "dens_mol_sol_comp": ConstantProperties.Constant,
            "visc_d_phase_comp": {"Sol": ConstantProperties.Constant},
            'valid_phase_types': PT.solidPhase,
            "parameter_data": {
                "mw": (66.37, pyunits.g / pyunits.mol),
                "cp_mol_sol_comp_coeff": (68.27, pyunits.J/pyunits.mol/pyunits.K), #Cp weighted average by composition of constituents
                "dens_mol_sol_comp_coeff": (2960.415544, pyunits.mol/pyunits.m**3), #ignore
                "enth_mol_form_sol_comp_ref": (0, pyunits.kJ/pyunits.mol),          #ignore
                "enrt_mol_form_sol_comp_ref": (158.1, pyunits.J/pyunits.mol/pyunits.K), #ignoer
                "visc_d_Sol_comp_coeff": (3.2833e-05, pyunits.Pa*pyunits.s)  
            },   

        },        
    },

# Specifying phases
"phases": {
    # "Vap": {"type": VaporPhase, "equation_of_state": Ideal},#Pv=nT
    "Sol": {"type": SolidPhase, "equation_of_state": Ideal},

},
# Set base units of measurement
"base_units": {
    "time": pyunits.s,
    "length": pyunits.m,
    "mass": pyunits.kg,
    "amount": pyunits.mol,
    "temperature": pyunits.K,
},
# Specifying state definition
"state_definition": FTPx,
"state_bounds": {
    "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
    "temperature": (273.15, 300, 2500, pyunits.K),
    "pressure": (5e3, 1e5, 1e6, pyunits.Pa),
},
"pressure_ref": (1e5, pyunits.Pa),
"temperature_ref": (300, pyunits.K),

"include_enthalpy_of_formation":(False)
}
