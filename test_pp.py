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
from idaes.models.properties.modular_properties.transport_properties.viscosity_wilke import ViscosityWilke
from idaes.models.properties.modular_properties.pure import ChapmanEnskog
from idaes.models.properties.modular_properties.pure.ChapmanEnskog import collision_integral_neufeld_callback 
from idaes.models.properties.modular_properties.pure import NIST
from idaes.core import PhaseType as PT
from idaes.models_extra.power_generation.properties.flue_gas_ideal import FlueGasStateBlock
from idaes.models.properties.modular_properties.transport_properties.viscosity_wilke import wilke_phi_ij_callback, herring_zimmer_phi_ij_callback

# Set up logger
_log = idaeslog.getLogger(__name__)

configuration = {
    # Specifying components
    "components": {
        
        "O2": {
            "type": Component,
            "elemental_composition": {"O":2},
            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            "visc_d_phase_comp": {"Vap": ChapmanEnskogLennardJones},
            'valid_phase_types': PT.vaporPhase,
            "parameter_data": {
                "mw": (31.9988, pyunits.g / pyunits.mol),  # [4]
                "pressure_crit": (50.43e5, pyunits.Pa),  # [8]
                "temperature_crit": (154.58, pyunits.K),  # [8]
                "cp_mol_ig_comp_coeff": { # valid range 100 K - 700 K
                    "A": (30.03235	, pyunits.J / pyunits.mol / pyunits.K),  # [4] #valid range  700-2000
                    "B": (8.772972, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1),
                    "C": (-3.988133, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2),
                    "D": (0.788313, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3),
                    "E": (-0.741599, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2),
                    "F": (-11.32468, pyunits.kJ / pyunits.mol),
                    "G": (236.1663, pyunits.J / pyunits.mol /pyunits.K),
                    "H": (0, pyunits.kJ / pyunits.mol)
                },
                "enth_mol_form_vap_comp_ref": (0, pyunits.kJ / pyunits.mol),  # [4]
                # "visc_d_Vap_comp_coeff": (3.8642e-05	, pyunits.Pa*pyunits.s)
                "lennard_jones_sigma": (2.641, pyunits.angstrom),
                "lennard_jones_epsilon_reduced": (809.1, pyunits.K),

            },
        },     
    },

# Specifying phases
"phases": {
    "Vap": {"type": VaporPhase, 
            "equation_of_state": Ideal,
            "visc_d_phase": ViscosityWilke,
            "transport_property_options": {
        "viscosity_phi_ij_callback": wilke_phi_ij_callback,
      }
            },#Pv=nT
    # "Sol": {"type": SolidPhase, "equation_of_state": Ideal},

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
