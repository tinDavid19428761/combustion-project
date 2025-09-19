#modular_biomass_rxn
from pyomo.environ import units as pyunits

from idaes.core import MaterialFlowBasis

from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn
from idaes.models.properties.modular_properties.reactions.rate_constant import arrhenius
from idaes.models.properties.modular_properties.reactions.rate_forms import power_law_rate

from idaes.models.properties.modular_properties.base.utility import ConcentrationForm

rxn_configuration = {
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "rate_reactions":{
        "Rcoal": {
            "stoichiometry": { #assumed 3wt% ash + ultimate analysis [C,H,O,N,S] from https://pmc.ncbi.nlm.nih.gov/articles/PMC9178619/
                            ("Vap", "coal"): -1,
                            ("Vap", "O2"): -0.8211,
                            ("Vap", "CO2"): 0.8082,
                            ("Vap", "H2O"): 0.01988,
                            ("Vap", "NO2"): 0.011134,
                            ("Vap", "SO2"): 0.036682,
                            ("Vap", "uncombustible"): 0.005905,
                            },
            "heat_of_reaction": constant_dh_rxn,
                "rate_constant": arrhenius,
                "rate_form": power_law_rate,
                "concentration_form": ConcentrationForm.moleFraction,
                "parameter_data": {
                    "dh_rxn_ref": (-376200.0747, pyunits.J/pyunits.mol),
                    "arrhenius_const": (1, pyunits.mol/pyunits.m**3/pyunits.s),
                    "energy_activation": (0, pyunits.J/pyunits.mol)},
        }
    }
}