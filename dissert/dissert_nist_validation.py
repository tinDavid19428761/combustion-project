import numpy as np

#stoich_reactor_test
from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    value,
    units as pyunits
)
#Todo add the four other unit operations
from idaes.models.unit_models import (Pump)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
# Import idaes logger to set output levels
import idaes.logger as idaeslog
from idaes.models.properties.modular_properties import GenericParameterBlock
# from  biomass_comb_pp import configuration 


# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import VaporPhase, Component
import idaes.logger as idaeslog

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ideal import Ideal

from idaes.models.properties.modular_properties.pure import NIST
from idaes.core import PhaseType as PT


# Set up logger
_log = idaeslog.getLogger(__name__)


#modular property package of just NIST ideal gases
configuration = {
    "include_enthalpy_of_formation":(False), 
    # Specifying components
    "components": {
        "CH4": {
            "type": Component,
            "elemental_composition": {"C": 1, "H": 4},
            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            'valid_phase_types': PT.vaporPhase,
            "parameter_data": {
                "mw": (16.0425, pyunits.g / pyunits.mol),  # [4]
                "pressure_crit": (46.1e5, pyunits.Pa),  # [[4]
                "temperature_crit": (190.6, pyunits.K),  # [4]
                "cp_mol_ig_comp_coeff": { #valid range 298 K - 1300 K
                    "A": (-0.703029	, pyunits.J / pyunits.mol / pyunits.K),  # [4]
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

        "O2": {
            "type": Component,
            "elemental_composition": {"O":2},
            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            'valid_phase_types': PT.vaporPhase,
            "parameter_data": {
                "mw": (31.9988, pyunits.g / pyunits.mol),  # [4]
                "pressure_crit": (50.43e5, pyunits.Pa),  # [8]
                "temperature_crit": (154.58, pyunits.K),  # [8]
                "cp_mol_ig_comp_coeff": { # valid range 100 K - 700 K
                    "A": (31.32234	, pyunits.J / pyunits.mol / pyunits.K),  # [4]
                    "B": (-20.23531, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1),
                    "C": (57.86644, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2),
                    "D": (-36.50624, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3),
                    "E": (-0.007374, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2),
                    "F": (-8.903471, pyunits.kJ / pyunits.mol),
                    "G": (246.7945, pyunits.J / pyunits.mol /pyunits.K),
                    "H": (0, pyunits.kJ / pyunits.mol)
                    #valid range 700 K - 2000 K
                    # "A": (30.03235	, pyunits.J / pyunits.mol / pyunits.K),  # [4] #valid range  700-2000
                    # "B": (8.772972, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1),
                    # "C": (-3.988133, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2),
                    # "D": (0.788313, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3),
                    # "E": (-0.741599, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2),
                    # "F": (-11.32468, pyunits.kJ / pyunits.mol),
                    # "G": (236.1663, pyunits.J / pyunits.mol /pyunits.K),
                    # "H": (0, pyunits.kJ / pyunits.mol)
                },
                "enth_mol_form_vap_comp_ref": (0, pyunits.kJ / pyunits.mol),

            },
        },


        "CO2": {
            "type": Component,
            "elemental_composition": {"C":1,"O":2, },
            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            'valid_phase_types': PT.vaporPhase,
            "parameter_data": {
                "mw": (44.0095, pyunits.g / pyunits.mol),  # [4]
                "pressure_crit": (73.80e5, pyunits.Pa),  # [[6]
                "temperature_crit": (304.18, pyunits.K),  # [6]
                "cp_mol_ig_comp_coeff": { #valid range 298 K - 1200 K
                    "A": (24.99735	, pyunits.J / pyunits.mol / pyunits.K),  # [4]
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
           'valid_phase_types': PT.vaporPhase,
           "parameter_data": {
               "mw": (18.0153e-3, pyunits.kg / pyunits.mol),
               "cp_mol_ig_comp_coeff": { #valid range 500 K- 1700 K
                   "A": (30.09200,pyunits.J / pyunits.mol / pyunits.K,),  # [4] 
                   "B": (6.832514,pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1,),
                   "C": (6.793435,pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2,),
                   "D": (-2.534480,pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3,),
                   "E": (0.082139,pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2,),
                   "F": (-250.8810, pyunits.kJ / pyunits.mol),
                   "G": (223.3967, pyunits.J / pyunits.mol / pyunits.K),
                   "H": (-241.8264, pyunits.kJ / pyunits.mol),
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
                "mw": (28.0134, pyunits.g / pyunits.mol),  # [4]
                "pressure_crit": (33.978e5, pyunits.Pa),  # [[7]
                "temperature_crit": (126.19, pyunits.K),  # [7]
                # "cp_mol_ig_comp_coeff": { #valid range 100 K to 500 K
                #     "A": (28.98641		, pyunits.J / pyunits.mol / pyunits.K),  # [4]
                #     "B": (1.853978, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1),
                #     "C": (-9.647459	, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2),
                #     "D": (16.63537	, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3),
                #     "E": (0.000117, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2),
                #     "F": (-8.671914, pyunits.kJ / pyunits.mol),
                #     "G": (226.4168, pyunits.J / pyunits.mol /pyunits.K),
                #     "H": (0, pyunits.kJ / pyunits.mol)
                # },
                "cp_mol_ig_comp_coeff": { #valid range 500 K to 2000 K
                    "A": (19.50583	, pyunits.J / pyunits.mol / pyunits.K),  # [4]
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
    "Vap": {"type": VaporPhase, "equation_of_state": Ideal},#Pv=nT

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


m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

m.fs.methane = GenericParameterBlock(**configuration)

m.fs.react = Pump(
    property_package = m.fs.methane,
)


m.fs.react.inlet.mole_frac_comp[0,"O2"].fix(1e-20)
m.fs.react.inlet.mole_frac_comp[0,"N2"].fix(1e-20)
m.fs.react.inlet.mole_frac_comp[0,"CO2"].fix(1e-20)
m.fs.react.inlet.mole_frac_comp[0,"H2O"].fix(1e-20)
m.fs.react.inlet.mole_frac_comp[0,"CH4"].fix(1e-20)
m.fs.react.inlet.temperature.fix(1000)
m.fs.react.inlet.pressure.fix(1e5)
m.fs.react.inlet.flow_mol.fix(10)
m.fs.react.deltaP.fix(30000)
m.fs.react.efficiency_pump.fix(0.8)


print(degrees_of_freedom(m))
assert degrees_of_freedom(m) == 0
m.fs.react.initialize(outlvl=idaeslog.INFO)
solver=SolverFactory("ipopt")
status=solver.solve(m,tee=True)

temps = list(range(300,2001,100))
temps.insert(0,298.15)
Cps = np.zeros((19,5)) #(row,col)
Cpstr = ""
compounds = ["O2","N2","CO2","H2O","CH4"]

for i in range(5):
    for j in range(19):
        m.fs.react.inlet.mole_frac_comp[0,compounds[i]].fix(1)
        m.fs.react.inlet.mole_frac_comp[0,compounds[i-1]].fix(1e-20)
        m.fs.react.inlet.temperature.fix(temps[j])
        solver=SolverFactory("ipopt")
        status=solver.solve(m,tee=True)
        Cps[j][i] = value(m.fs.react.control_volume.properties_in[0].cp_mol)


for j in range(19):
    for i in range(5):
        Cpstr+=str("{:.3f}".format(Cps[j][i]))
        Cpstr+=","
    Cpstr+="\n"

print(Cps)
print(Cpstr)