"""
Benzene-Toluene phase equilibrium package using ideal liquid and vapor.

Example property package using the Generic Property Package Framework.
This exmample shows how to set up a property package to do benzene-toluene
phase equilibrium in the generic framework using ideal liquid and vapor
assumptions along with methods drawn from the pre-built IDAES property
libraries.
"""
# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component, SolidPhase
import idaes.logger as idaeslog

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.pure import ConstantProperties
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


# ---------------------------------------------------------------------
# Configuration dictionary for an ideal Benzene-Toluene system

# Data Sources:
# [1] The Properties of Gases and Liquids (1987)
#     4th edition, Chemical Engineering Series - Robert C. Reid
# [2] Perry's Chemical Engineers' Handbook 7th Ed.
# [3] Engineering Toolbox, https://www.engineeringtoolbox.com
#     Retrieved 1st December, 2019
#[4]NIST
#[5] wagner, Ewers, et al.,1976
#[6] Suehiro, Nakajima, et al., 1996
#[7] Jacobsen, Stewart, et al., 1986
# [8] 	Cardoso, 1915

configuration = {
    # Specifying components
    "components": {
        # https://github.com/IDAES/idaes-pse/blob/main/idaes/models/properties/modular_properties/pure/ConstantProperties.py
        #constant[solid]properties (pure)
        # parameters: : cp_mol, enth_mol_form, entr_mol_form, dens_mol()
        "carbon": { 
            "type": Component,
            "elemental_composition": {"C":1},
            "enth_mol_sol_comp": ConstantProperties.Constant,
            "dens_mol_sol_comp": ConstantProperties.Constant,
            'valid_phase_types': PT.solidPhase,

            "parameter_data": {
                "cp_mol_sol_comp_coeff": (10, pyunits.J/pyunits.mol/pyunits.K),
                # 1100-1800 kg/m3 to mol/m3 = density*1000/12 ~=150,000 mol/m3, (for coal)
                "dens_mol_sol_comp_coeff": (150000, pyunits.mol/pyunits.m**3),
                "enth_mol_form_sol_comp_ref": (0, pyunits.kJ/pyunits.mol),
                "enrt_mol_form_sol_comp_ref": (158.1, pyunits.J/pyunits.mol/pyunits.K)
            },      
        },

        "CH4": {
            "type": Component,
            "elemental_composition": {"C": 1, "H": 4},
            "enth_mol_ig_comp": NIST,
            'valid_phase_types': PT.vaporPhase,
            "parameter_data": {
                "mw": (16.0425, pyunits.g / pyunits.mol),  # [4]
                "pressure_crit": (46.1e5, pyunits.Pa),  # [[4]
                "temperature_crit": (190.6, pyunits.K),  # [4]
                "cp_mol_ig_comp_coeff": {
                    "A": (-0.703029	, pyunits.J / pyunits.mol / pyunits.K),  # [4]
                    "B": (108.4773, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1),
                    "C": (-42.52157, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2),
                    "D": (5.862788, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3),
                    "E": (0.678565, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2),
                    "F": (-76.84376, pyunits.kJ / pyunits.mol),
                    "G": (158.7163, pyunits.J / pyunits.mol /pyunits.K),
                    "H": (-74.87310, pyunits.kJ / pyunits.mol)
                },
                

                "enth_mol_form_vap_comp_ref": (-74.87310, pyunits.kJ / pyunits.mol),  # [4]
            },
        },
        "O2": {
            "type": Component,
            "elemental_composition": {"O":2},
            "enth_mol_ig_comp": NIST,
            'valid_phase_types': PT.vaporPhase,
            "parameter_data": {
                "mw": (31.9988, pyunits.g / pyunits.mol),  # [4]
                "pressure_crit": (50.43e5, pyunits.Pa),  # [8]
                "temperature_crit": (154.58, pyunits.K),  # [8]
                "cp_mol_ig_comp_coeff": {
                    "A": (31.32234	, pyunits.J / pyunits.mol / pyunits.K),  # [4]
                    "B": (-20.23531, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1),
                    "C": (57.86644, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2),
                    "D": (-36.50624, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3),
                    "E": (-0.007374, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2),
                    "F": (-8.903471, pyunits.kJ / pyunits.mol),
                    "G": (246.7945, pyunits.J / pyunits.mol /pyunits.K),
                    "H": (0, pyunits.kJ / pyunits.mol)
                },
                
                "enth_mol_form_vap_comp_ref": (0, pyunits.kJ / pyunits.mol),  # [4]
            },
        },
        "CO2": {
            "type": Component,
            "elemental_composition": {"O":2, "C":1},
            "enth_mol_ig_comp": NIST,
            'valid_phase_types': PT.vaporPhase,
            "parameter_data": {
                "mw": (44.0095, pyunits.g / pyunits.mol),  # [4]
                "pressure_crit": (73.80e5, pyunits.Pa),  # [[6]
                "temperature_crit": (304.18, pyunits.K),  # [6]
                "cp_mol_ig_comp_coeff": {
                    "A": (24.99735	, pyunits.J / pyunits.mol / pyunits.K),  # [4]
                    "B": (55.18696, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1),
                    "C": (-33.69137, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2),
                    "D": (7.948387, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3),
                    "E": (-0.136638, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2),
                    "F": (-403.6075, pyunits.kJ / pyunits.mol),
                    "G": (228.2431, pyunits.J / pyunits.mol /pyunits.K),
                    "H": (-393.5224, pyunits.kJ / pyunits.mol)
                },
                
                "enth_mol_form_vap_comp_ref": (-393.5224, pyunits.kJ / pyunits.mol),  # [4]
            },
        },
        "CO": {
            "type": Component,
            "elemental_composition": {"O":1, "C":1},
            "enth_mol_ig_comp": NIST,
            'valid_phase_types': PT.vaporPhase,
            "parameter_data": {
                "mw": (28.0101, pyunits.g / pyunits.mol),  # [4]
                "pressure_crit": (34.9875e5, pyunits.Pa),  # [[6]
                "temperature_crit": (134.45, pyunits.K),  # [6]
                "cp_mol_ig_comp_coeff": { #(valid range: 298-1300K)
                    "A": (25.56759	, pyunits.J / pyunits.mol / pyunits.K),  # [4]
                    "B": (6.096130, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1),
                    "C": (4.054656, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2),
                    "D": (-2.671301	, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3),
                    "E": (0.131021, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2),
                    "F": (-118.0089, pyunits.kJ / pyunits.mol),
                    "G": (227.3665, pyunits.J / pyunits.mol /pyunits.K),
                    "H": (-110.5271, pyunits.kJ / pyunits.mol)
                },
                
                "enth_mol_form_vap_comp_ref": (-110.5271, pyunits.kJ / pyunits.mol),  # [4]
            },
        },
        "H2O": {
           "type": Component,
           "elemental_composition":{"H":2,"O":1},
           "enth_mol_ig_comp": NIST,
           "pressure_sat_comp": NIST,
           'valid_phase_types': PT.vaporPhase,
           "parameter_data": {
               "mw": (18.0153e-3, pyunits.kg / pyunits.mol),  # [4]
               "pressure_crit": (220.64e5, pyunits.Pa),  # [4]
               "temperature_crit": (647, pyunits.K),  # [4]
               "cp_mol_ig_comp_coeff": {
                   "A": (
                       30.09200,
                       pyunits.J / pyunits.mol / pyunits.K,
                   ),  # [4] temperature range 500 K- 1700 K
                   "B": (
                       6.832514,
                       pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1,
                   ),
                   "C": (
                       6.793435,
                       pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2,
                   ),
                   "D": (
                       -2.534480,
                       pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3,
                   ),
                   "E": (
                       0.082139,
                       pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2,
                   ),
                   "F": (-250.8810, pyunits.kJ / pyunits.mol),
                   "G": (223.3967, pyunits.J / pyunits.mol / pyunits.K),
                   "H": (-241.83, pyunits.kJ / pyunits.mol),
               },
               "enth_mol_form_vap_comp_ref": (-241830, pyunits.J / pyunits.mol),  # [4]
               "pressure_sat_comp_coeff": {
                   "A": (4.6543, None),  # [4], temperature range 255.9 K - 373 K
                   "B": (1435.264, pyunits.K),
                   "C": (-64.848, pyunits.K),
               },
           },
       },
        "N2": {
            "type": Component,
            "elemental_composition": {"N":2},
            "enth_mol_ig_comp": NIST,
            'valid_phase_types': PT.vaporPhase,
            "parameter_data": {
                "mw": (28.0134, pyunits.g / pyunits.mol),  # [4]
                "pressure_crit": (33.978e5, pyunits.Pa),  # [[7]
                "temperature_crit": (126.19, pyunits.K),  # [7]
                "cp_mol_ig_comp_coeff": {
                    "A": (28.98641		, pyunits.J / pyunits.mol / pyunits.K),  # [4]
                    "B": (1.853978, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1),
                    "C": (-9.647459	, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2),
                    "D": (16.63537	, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3),
                    "E": (0.000117, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2),
                    "F": (-8.671914, pyunits.kJ / pyunits.mol),
                    "G": (226.4168, pyunits.J / pyunits.mol /pyunits.K),
                    "H": (0, pyunits.kJ / pyunits.mol)
                },
                "enth_mol_form_vap_comp_ref": (0, pyunits.kJ / pyunits.mol),  # [4]
            },
        },
    },

    # Specifying phases
"phases": {
    "Vap": {"type": VaporPhase, "equation_of_state": Ideal},#Pv=nT
    "Sol": {"type": SolidPhase, "equation_of_state": Ideal}

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

}
