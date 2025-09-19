#modular first order reaction
from pyomo.environ import units as pyunits

from idaes.core import MaterialFlowBasis

from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn
from idaes.models.properties.modular_properties.reactions.rate_constant import arrhenius
from idaes.models.properties.modular_properties.reactions.rate_forms import power_law_rate

from idaes.models.properties.modular_properties.base.utility import ConcentrationForm

config_dict = {
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "rate_reactions": {
        "R1": {"stoichiometry": {("Liq", "H2O2"): -2,
                                 ("Vap", "H2O"): 2,
                                 ("Vap", "O2"):  1},
               "heat_of_reaction": constant_dh_rxn,
               "rate_constant": arrhenius,
               "rate_form": power_law_rate,
               "concentration_form": ConcentrationForm.moleFraction,
               "parameter_data": {
                   "dh_rxn_ref": (-0, pyunits.J/pyunits.mol),#zero, calc by enth_forms instead
                   "arrhenius_const": (17.70574, pyunits.mol/pyunits.m**3/pyunits.s),
                   "energy_activation": (30000, pyunits.J/pyunits.mol)}}},
}