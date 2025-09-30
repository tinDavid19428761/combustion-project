"""
Property package for the combustion of biomass in air
"""
# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component, SolidPhase
import idaes.logger as idaeslog

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.pure import ConstantProperties, ChapmanEnskog
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    IdealBubbleDew,
)
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity
from idaes.models.properties.modular_properties.pure import Perrys
from idaes.models.properties.modular_properties.pure import RPP4
from idaes.models.properties.modular_properties.pure import NIST
from idaes.core import PhaseType as PT


# Set up logger
_log = idaeslog.getLogger(__name__)

configuration = {
    # Specifying components
    "components": {
        "biomass": { #woody biomass
            "type": Component,
            "elemental_composition": {"C":6, "H": 10, "O": 5}, #cellulose composition C6-H10-O5
            "enth_mol_sol_comp": ConstantProperties.Constant,
            "cp_mol_sol_comp": ConstantProperties.Constant,
            "dens_mol_sol_comp": ConstantProperties.Constant,
            "visc_d_phase_comp": {"Sol": ConstantProperties.Constant},
            'valid_phase_types': PT.solidPhase,
            "parameter_data": {
                "mw": (162.1394, pyunits.g / pyunits.mol),
                "cp_mol_sol_comp_coeff": (243.2091, pyunits.J/pyunits.mol/pyunits.K),
                "dens_mol_sol_comp_coeff": (2960.415544, pyunits.mol/pyunits.m**3), 
                "enth_mol_form_sol_comp_ref": (0, pyunits.kJ/pyunits.mol),
                "enrt_mol_form_sol_comp_ref": (158.1, pyunits.J/pyunits.mol/pyunits.K),
                "visc_d_Sol_comp_coeff": (3.2833e-05, pyunits.Pa*pyunits.s)  
            }, 
        },
        "uncombustible": { #wood ash
            "type": Component,
            "elemental_composition": {"uncombustible":1}, 
            "cp_mol_sol_comp": ConstantProperties.Constant,
            "enth_mol_sol_comp": ConstantProperties.Constant,
            "dens_mol_sol_comp": ConstantProperties.Constant,
            'valid_phase_types': PT.solidPhase,
            "parameter_data": {
                "mw": (66.37, pyunits.g / pyunits.mol),
                "cp_mol_sol_comp_coeff": (68.27, pyunits.J/pyunits.mol/pyunits.K), 
                "dens_mol_sol_comp_coeff": (2960.415544, pyunits.mol/pyunits.m**3), 
                "enth_mol_form_sol_comp_ref": (0, pyunits.kJ/pyunits.mol),       
                "enrt_mol_form_sol_comp_ref": (158.1, pyunits.J/pyunits.mol/pyunits.K), 
                "visc_d_Sol_comp_coeff": (3.2833e-05, pyunits.Pa*pyunits.s)  
            },   
        },
        "O2": {
            "type": Component,
            "elemental_composition": {"O":2},
            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            'valid_phase_types': PT.vaporPhase,
            "parameter_data": {
                "mw": (31.9988, pyunits.g / pyunits.mol), 
                "cp_mol_ig_comp_coeff": {  #valid range  700-2000
                    "A": (30.03235	, pyunits.J / pyunits.mol / pyunits.K), 
                    "B": (8.772972, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1),
                    "C": (-3.988133, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2),
                    "D": (0.788313, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3),
                    "E": (-0.741599, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2),
                    "F": (-11.32468, pyunits.kJ / pyunits.mol),
                    "G": (236.1663, pyunits.J / pyunits.mol /pyunits.K),
                    "H": (0, pyunits.kJ / pyunits.mol)
                },
            },
        },
        "CH4": {
            "type": Component,
            "elemental_composition": {"C": 1, "H": 4},
            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            'valid_phase_types': PT.vaporPhase,
            "parameter_data": {
                "mw": (16.0425, pyunits.g / pyunits.mol),  
                "cp_mol_ig_comp_coeff": { #valid range 298 K - 1300 K
                    "A": (-0.703029	, pyunits.J / pyunits.mol / pyunits.K),
                    "B": (108.4773, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1),
                    "C": (-42.52157, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2),
                    "D": (5.862788, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3),
                    "E": (0.678565, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2),
                    "F": (-76.84376, pyunits.kJ / pyunits.mol),
                    "G": (158.7163, pyunits.J / pyunits.mol /pyunits.K),
                    "H": (-74.87310, pyunits.kJ / pyunits.mol),
                },
            },
        },
        "CO2": {
            "type": Component,
            "elemental_composition": {"C":1,"O":2, },
            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            'valid_phase_types': PT.vaporPhase,
            "parameter_data": {
                "mw": (44.0095, pyunits.g / pyunits.mol), 
                "pressure_crit": (73.80e5, pyunits.Pa),  
                "temperature_crit": (304.18, pyunits.K),  
                "cp_mol_ig_comp_coeff": { #valid range 298 K - 1200 K
                    "A": (24.99735	, pyunits.J / pyunits.mol / pyunits.K), 
                    "B": (55.18696, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1),
                    "C": (-33.69137, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2),
                    "D": (7.948387, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3),
                    "E": (-0.136638, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2),
                    "F": (-403.6075, pyunits.kJ / pyunits.mol),
                    "G": (228.2431, pyunits.J / pyunits.mol /pyunits.K),
                    "H": (-393.5224, pyunits.kJ / pyunits.mol),
                },
            },
        },
        "H2O": {
           "type": Component,
           "elemental_composition":{"H":2,"O":1},
           "enth_mol_ig_comp": NIST,
           "cp_mol_ig_comp": NIST,
           "pressure_sat_comp": NIST,
           'valid_phase_types': PT.vaporPhase,
           "parameter_data": {
               "mw": (18.0153e-3, pyunits.kg / pyunits.mol), 
               "cp_mol_ig_comp_coeff": { #valid range 500 K- 1700 K
                   "A": (30.09200,pyunits.J / pyunits.mol / pyunits.K,),  
                   "B": (6.832514,pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1,),
                   "C": (6.793435,pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2,),
                   "D": (-2.534480,pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3,),
                   "E": (0.082139,pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2,),
                   "F": (-250.8810, pyunits.kJ / pyunits.mol),
                   "G": (223.3967, pyunits.J / pyunits.mol / pyunits.K),
                   "H": (-241.8264, pyunits.kJ / pyunits.mol),
               },
               "pressure_sat_comp_coeff": {
                   "A": (4.6543, None), #temperature range 255.9 K - 373 K
                   "B": (1435.264, pyunits.K),
                   "C": (-64.848, pyunits.K),
               },
           },
       },
        "N2": {
            "type": Component,
            "elemental_composition": {"N":2},
            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            'valid_phase_types': PT.vaporPhase,
            "parameter_data": {
                "mw": (28.0134, pyunits.g / pyunits.mol),  
                "cp_mol_ig_comp_coeff": { #valid range 500 K to 2000 K
                    "A": (19.50583	, pyunits.J / pyunits.mol / pyunits.K), 
                    "B": (19.88705, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1),
                    "C": (-8.598535	, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2),
                    "D": (1.369784	, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3),
                    "E": (0.527601, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2),
                    "F": (-4.935202, pyunits.kJ / pyunits.mol),
                    "G": (212.3900, pyunits.J / pyunits.mol /pyunits.K),
                    "H": (0, pyunits.kJ / pyunits.mol)
                },
            },
        },
    },

# Specifying phases
"phases": {
    "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
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

"include_enthalpy_of_formation":(False), 
}
