"""
reaction package for the combustion of coal in air
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
                                              'coal',
                                              'uncombustible',
                                              'CH4'
                                              ])

        # Reaction Index
        self.rate_reaction_idx = Set(initialize=["Rcoal","RCH4"])
        self.uncombs_set = Set(initialize=["Rcoal",])
    
        self.reaction_set = Set(initialize=[("Rcoal", "Vap", "H2O"),
                                            ("Rcoal", "Vap", "CO2"),
                                            ("Rcoal", "Vap", "O2"),
                                            ("Rcoal", "Sol", "coal"),
                                            ("Rcoal", "Vap", "N2"),
                                            ("Rcoal", "Sol", "uncombustible"),
                                            ("Rcoal", "Vap", "CH4"),

                                            ("RCH4", "Vap", "H2O"),
                                            ("RCH4", "Vap", "CO2"),
                                            ("RCH4", "Vap", "O2"),
                                            ("RCH4", "Sol", "coal"),
                                            ("RCH4", "Vap", "N2"),
                                            ("RCH4", "Sol", "uncombustible"),
                                            ("RCH4", "Vap", "CH4"),
                                            ])

        # Reaction Stoichiometry
        self.rate_reaction_stoichiometry = Var(self.reaction_set, initialize={
                                            ("Rcoal", "Vap", "H2O"): 0.161580373,
                                            ("Rcoal", "Vap", "CO2"): 0.715621569,
                                            ("Rcoal", "Vap", "O2"): -0.643121569,
                                            ("Rcoal", "Sol", "coal"): -1,
                                            ("Rcoal", "Vap", "N2"): 0.005,
                                            ("Rcoal", "Sol", "uncombustible"): 0.03,
                                            ("Rcoal", "Vap", "CH4"): 0,

                                            ("RCH4", "Vap", "H2O"): 2,
                                            ("RCH4", "Vap", "CO2"): 1,
                                            ("RCH4", "Vap", "O2"): -2,
                                            ("RCH4", "Sol", "coal"): 0,
                                            ("RCH4", "Vap", "N2"): 0,
                                            ("RCH4", "Sol", "uncombustible"): 0,
                                            ("RCH4", "Vap", "CH4"): -1,
                                            })
        self.rate_reaction_stoichiometry.fix()
        
        # self.reactant_list=Set(initialize=["coal","O2",'CH4'])
        self.limit_reactant_dict = Param(self.rate_reaction_idx, initialize={
            "Rcoal": "coal",
            "RCH4": "CH4",
        },
        within=Any)

        dh_rxn_dict = {"Rcoal": -328066.056, 
                       "RCH4": -802125 
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
    

