#H2O2 CSTr test

import numpy as np

#stoich_reactor_test
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
# from pytest import approx

#Todo add the four other unit operations
from idaes.models.unit_models import (
Mixer,
StoichiometricReactor,
CSTR,
Heater,
HeatExchanger,
Separator,
Flash,
Heater,
Pump
)
from idaes.models_extra.power_generation.unit_models.helm.phase_separator import HelmPhaseSeparator

from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
# Import idaes logger to set output levels
import idaes.logger as idaeslog
from idaes.models.properties.modular_properties import GenericParameterBlock, GenericReactionParameterBlock
# from  biomass_comb_pp import configuration 


#helmholtz import for water
from idaes.models.properties.general_helmholtz import (
        HelmholtzParameterBlock,
        HelmholtzThermoExpressions,
        AmountBasis,
        PhaseType,
        StateVars
    )
from idaes.models.unit_models.heat_exchanger import delta_temperature_amtd_callback, HX0DInitializer
from idaes.models.unit_models.separator import SplittingType, EnergySplittingType
from idaes.core.util.model_diagnostics import (
    DiagnosticsToolbox,
)

# from custom_combustion_reactor import MultiCombReactor
import unittest


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
    "include_enthalpy_of_formation":(False), #put this at the top
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
            "visc_d_phase_comp": {"Vap": ConstantProperties.Constant},
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
                    # "A": (30.03235	, pyunits.J / pyunits.mol / pyunits.K),  # [4] #valid range  700-2000
                    # "B": (8.772972, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1),
                    # "C": (-3.988133, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2),
                    # "D": (0.788313, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3),
                    # "E": (-0.741599, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2),
                    # "F": (-11.32468, pyunits.kJ / pyunits.mol),
                    # "G": (236.1663, pyunits.J / pyunits.mol /pyunits.K),
                    # "H": (0, pyunits.kJ / pyunits.mol)
                },
                "enth_mol_form_vap_comp_ref": (0, pyunits.kJ / pyunits.mol),  # [4]
                "visc_d_Vap_comp_coeff": (3.8642e-05	, pyunits.Pa*pyunits.s)

            },
        },


        "CO2": {
            "type": Component,
            "elemental_composition": {"C":1,"O":2, },
            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            "visc_d_phase_comp": {"Vap": ConstantProperties.Constant},
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
                "enth_mol_form_vap_comp_ref": (0, pyunits.kJ / pyunits.mol),  # [4]
                "visc_d_Vap_comp_coeff": (3.2115e-05	, pyunits.Pa*pyunits.s)

            },
        },

        # "CO": {
        #     "type": Component,
        #     "elemental_composition": {"O":1, "C":1},
        #     "enth_mol_ig_comp": NIST,
        #     "cp_mol_ig_comp": NIST,
        #     "visc_d_phase_comp": {"Vap": ConstantProperties.Constant},
        #     'valid_phase_types': PT.vaporPhase,
        #     "parameter_data": {
        #         "mw": (28.0101, pyunits.g / pyunits.mol),  # [4]
        #         "pressure_crit": (34.9875e5, pyunits.Pa),  # [[6]
        #         "temperature_crit": (134.45, pyunits.K),  # [6]
        #         "cp_mol_ig_comp_coeff": { #(valid range: 298-1300K)
        #             "A": (25.56759	, pyunits.J / pyunits.mol / pyunits.K),  # [4] 
        #             "B": (6.096130, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1),
        #             "C": (4.054656, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2),
        #             "D": (-2.671301	, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3),
        #             "E": (0.131021, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2),
        #             "F": (-118.0089, pyunits.kJ / pyunits.mol),
        #             "G": (227.3665, pyunits.J / pyunits.mol /pyunits.K),
        #             "H": (-110.5271, pyunits.kJ / pyunits.mol),
        #         },
        #         "enth_mol_form_vap_comp_ref": (0, pyunits.kJ / pyunits.mol),  # [4]
        #         "visc_d_Vap_comp_coeff": (2.5942e-05	, pyunits.Pa*pyunits.s)

        #     },
        # },
        "H2O": {
           "type": Component,
           "elemental_composition":{"H":2,"O":1},
           "enth_mol_ig_comp": NIST,
           "cp_mol_ig_comp": NIST,
           "pressure_sat_comp": NIST,
           "visc_d_phase_comp": {"Vap": ConstantProperties.Constant},
           'valid_phase_types': PT.vaporPhase,
           "parameter_data": {
               "mw": (18.0153e-3, pyunits.kg / pyunits.mol),  # [4]
               "pressure_crit": (220.64e5, pyunits.Pa),  # [4]
               "temperature_crit": (647, pyunits.K),  # [4]
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
               "enth_mol_form_vap_comp_ref": (0, pyunits.J / pyunits.mol),  # [4]

               "pressure_sat_comp_coeff": {
                   "A": (4.6543, None),  # [4], temperature range 255.9 K - 373 K
                   "B": (1435.264, pyunits.K),
                   "C": (-64.848, pyunits.K),
               },
               "visc_d_Vap_comp_coeff": (2.8564e-05	, pyunits.Pa*pyunits.s)
           },
       },

        "N2": {
            "type": Component,
            "elemental_composition": {"N":2},
            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            "visc_d_phase_comp": {"Vap": ConstantProperties.Constant},
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
                "enth_mol_form_vap_comp_ref": (0, pyunits.kJ / pyunits.mol),  # [4]
                "visc_d_Vap_comp_coeff": (3.2833e-05, pyunits.Pa*pyunits.s)
            },
        },
        "NO": {
            "type": Component,
            "elemental_composition": {"N":1,"O":2, },
            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            'valid_phase_types': PT.vaporPhase,
            "parameter_data": {
                "mw": (46.0055, pyunits.g / pyunits.mol),  # [4]
                "pressure_crit": (73.80e5, pyunits.Pa),  # [[6]
                "temperature_crit": (304.18, pyunits.K),  # [6]
                "cp_mol_ig_comp_coeff": { #valid range 298 K - 1200 K
                    "A": (23.83491	, pyunits.J / pyunits.mol / pyunits.K),  # [4]
                    "B": (12.58878, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1),
                    "C": (-1.139011, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2),
                    "D": (-1.497459, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3),
                    "E": (0.214194, pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2),
                    "F": (83.35783, pyunits.kJ / pyunits.mol),
                    "G": (237.1219, pyunits.J / pyunits.mol /pyunits.K),
                    "H": (90.29114, pyunits.kJ / pyunits.mol),
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
m.fs.react.inlet.mole_frac_comp[0,"NO"].fix(1e-20)
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
Cps = np.zeros((19,6)) #(row,col)
Cpstr = ""
compounds = ["O2","N2","CO2","H2O","NO","CH4"]

for i in range(6):
    for j in range(19):
        m.fs.react.inlet.mole_frac_comp[0,compounds[i]].fix(1)
        m.fs.react.inlet.mole_frac_comp[0,compounds[i-1]].fix(1e-20)
        m.fs.react.inlet.temperature.fix(temps[j])
        solver=SolverFactory("ipopt")
        status=solver.solve(m,tee=True)
        Cps[j][i] = value(m.fs.react.control_volume.properties_in[0].cp_mol)
    # m.fs.react.initialize(outlvl=idaeslog.INFO)

for j in range(19):
    for i in range(6):
        Cpstr+=str("{:.3f}".format(Cps[j][i]))
        Cpstr+=","
    Cpstr+="\n"

print(Cps)
print(Cpstr)
# m.fs.react.report()
# print(value(m.fs.react.control_volume.properties_in[0].cp_mol))
# print(value(m.fs.react.control_volume.properties_in[0].enth_mol))
# m.fs.react.inlet.controlVolume0d

# Cpstr += str(value(m.fs.react.control_volume.properties_in[0].cp_mol))