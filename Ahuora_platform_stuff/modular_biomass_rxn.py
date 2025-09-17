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
        "Rbiomass": {
            "stoichiometry": {("Sol", "biomass"): -1,
                            ("Vap", "O2"): -6,
                            ("Vap", "CO2"): 6,
                            ("Vap", "H2O"): 5,
                            # ("Sol", "uncombustible"): 0.01,
                            },
            "heat_of_reaction": constant_dh_rxn,
                "rate_constant": arrhenius,
                "rate_form": power_law_rate,
                "concentration_form": ConcentrationForm.moleFraction,
                "parameter_data": {
                    "dh_rxn_ref": (-2749556.40, pyunits.J/pyunits.mol),
                    "arrhenius_const": (1, pyunits.mol/pyunits.m**3/pyunits.s),
                    "energy_activation": (0, pyunits.J/pyunits.mol)},
        }
    }
}