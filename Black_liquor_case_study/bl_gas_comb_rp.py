"""
reaction package for the combustion of biomass in air
"""
from pyomo.environ import Expression

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import (Constraint,
                           exp,
                           Set,
                           Var,
                           Param,
                           units as pyunits,
                           Reals,
                           Any)

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        ReactionParameterBlock,
                        ReactionBlockDataBase,
                        ReactionBlockBase)
from idaes.core.util.constants import Constants as const
from idaes.core.util.misc import add_object_reference

# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("MultiCombReactionParameterBlock")
class MultiCombReactionParameterData(ReactionParameterBlock):
    """
    Property Parameter Block Class
    Contains parameters and indexing sets associated with properties for
    superheated steam.
    """

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(MultiCombReactionParameterData, self).build()

        self._reaction_block_class = BMReactionBlock

        # List of valid phases in property package
        self.phase_list = Set(initialize=['Vap', 'Sol'])

        # Component list - a list of component identifiers
        self.component_list = Set(initialize=['H2O',
                                              'CO2',
                                              'O2',
                                              'N2',
                                              'BL',
                                              'uncombustible',
                                              'CH4'
                                              ])

        # Reaction Index
        self.rate_reaction_idx = Set(initialize=["Rbl","RCH4"])
        self.uncombs_set = Set(initialize=["Rbl",])
    
        self.reaction_set = Set(initialize=[("Rbl", "Vap", "H2O"),
                                            ("Rbl", "Vap", "CO2"),
                                            ("Rbl", "Vap", "O2"),
                                            ("Rbl", "Liq", "BL"),
                                            ("Rbl", "Vap", "N2"),
                                            ("Rbl", "Sol", "uncombustible"),
                                            ("Rbl", "Vap", "CH4"),

                                            ("RCH4", "Vap", "H2O"),
                                            ("RCH4", "Vap", "CO2"),
                                            ("RCH4", "Vap", "O2"),
                                            ("RCH4", "Liq", "BL"),
                                            ("RCH4", "Vap", "N2"),
                                            ("RCH4", "Sol", "uncombustible"),
                                            ("RCH4", "Vap", "CH4"),
                                            ])

        # Reaction Stoichiometry
        self.rate_reaction_stoichiometry = Var(self.reaction_set, initialize={
                                            ("Rbl", "Vap", "H2O"): 1-0.7-(0.875*0.3), # carryover water subtract(H2O consumed based on Co2 emit.)
                                            ("Rbl", "Vap", "CO2"): 0.3, #based on assumed black liquor emissions factor of 95.3 kgCO2/GJ   https://naturvardsverket.diva-portal.org/smash/get/diva2:1546963/FULLTEXT01.pdf      
                                            ("Rbl", "Vap", "O2"): -0.3+(0.15*0.875), #0.3 based on Co2 emit. add(Oxygen consumed/supplied for H2O)
                                            ("Rbl", "Liq", "BL"): -1,   
                                            ("Rbl", "Vap", "N2"): 0,
                                            # ("Rbl", "Sol", "uncombustible"): 0.7-0.3, #assume solids are retained
                                            ("Rbl", "Sol", "uncombustible"): 0.143745, #goal seek mass balance with other stoichs
                                            ("Rbl", "Vap", "CH4"): 0,

                                            ("RCH4", "Vap", "H2O"): 2,
                                            ("RCH4", "Vap", "CO2"): 1,
                                            ("RCH4", "Vap", "O2"): -2,
                                            ("RCH4", "Liq", "BL"): 0,
                                            ("RCH4", "Vap", "N2"): 0,
                                            ("RCH4", "Sol", "uncombustible"): 0,
                                            ("RCH4", "Vap", "CH4"): -1,
                                            })
        self.rate_reaction_stoichiometry.fix()
        
        # self.reactant_list=Set(initialize=["biomass","O2",'CH4'])
        self.limit_reactant_dict = Param(self.rate_reaction_idx, initialize={
            "Rbl": "BL",
            "RCH4": "CH4",
        },
        within=Any)

        dh_rxn_dict = {"Rbl": -135150, #
                       "RCH4": -802125 # (J/mol) engineering toolbox methane LHV
                       } 
        
        self.dh_rxn = Var(self.rate_reaction_idx, 
                          initialize = dh_rxn_dict,
                          domain=Reals,
                          doc="Heat of reaction",
                          units=pyunits.J/pyunits.mol)
        self.dh_rxn.fix()

        
        

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units({'time': pyunits.s,
                               'length': pyunits.m,
                               'mass': pyunits.kg,
                               'amount': pyunits.mol,
                               'temperature': pyunits.K})


class ReactionBlock(ReactionBlockBase):
    """
    This Class contains methods which should be applied to Reaction Blocks as a
    whole, rather than individual elements of indexed Reaction Blocks.
    """
    def initialize(blk, outlvl=0, **kwargs):
        '''
        Initialization routine for reaction package.
        Keyword Arguments:
            outlvl : sets output level of initialization routine
                     * 0 = no output (default)
                     * 1 = report after each step
        Returns:
            None
        '''
        if outlvl > 0:
            _log.info('{} Initialization Complete.'.format(blk.name))


@declare_process_block_class("BMReactionBlock", block_class=ReactionBlock)
class BMReactionBlockData(ReactionBlockDataBase):
    def build(self):
        """
        Callable method for Block construction
        """
        super(BMReactionBlockData, self).build()

        # Heat of reaction - no _ref as this is the actual property
        add_object_reference(
                self,
                "dh_rxn",
                self.config.parameters.dh_rxn) 

    def get_reaction_rate_basis(b):
        return MaterialFlowBasis.molar
    

