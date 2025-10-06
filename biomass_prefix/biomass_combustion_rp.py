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
                           Reals,
                           Any,
                           units as pyunits)

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


@declare_process_block_class("BMCombReactionParameterBlock")
class BMCombReactionParameterData(ReactionParameterBlock):
    """
    Property Parameter Block Class
    Contains parameters and indexing sets associated with properties for
    superheated steam.
    """

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(BMCombReactionParameterData, self).build()

        self._reaction_block_class = BMReactionBlock

        # List of valid phases in property package
        self.phase_list = Set(initialize=["Vap", "Sol"])

        # Component list - a list of component identifiers
        self.component_list = Set(initialize=["H2O",
                                              "CO2",
                                              "O2",
                                              "N2",
                                              "biomass",
                                              "ash",
                                              "coal",
                                              ])

        # Reaction Index
        self.rate_reaction_idx = Set(initialize=["R1","Rcoal","RCH4"])
        self.uncombs_set = Set(initialize=["R1","Rcoal"])
        


        self.reaction_set = Set(initialize=[("R1", "Vap", "H2O"),
                                            ("R1", "Vap", "CO2"),
                                            ("R1", "Vap", "O2"),
                                            ("R1", "Sol", "coal"),
                                            ("R1", "Sol", "biomass"),
                                            ("R1", "Vap", "N2"),
                                            ("R1", "Sol", "ash"),

                                            ("Rcoal", "Vap", "H2O"),
                                            ("Rcoal", "Vap", "CO2"),
                                            ("Rcoal", "Vap", "O2"),
                                            ("Rcoal", "Sol", "coal"),
                                            ("Rcoal", "Sol", "biomass"),
                                            ("Rcoal", "Vap", "N2"),
                                            ("Rcoal", "Sol", "ash"),

                                            ("RCH4", "Vap", "H2O"),
                                            ("RCH4", "Vap", "CO2"),
                                            ("RCH4", "Vap", "O2"),
                                            ("RCH4", "Vap", "CH4"),
                                            ])
        
        # default values for default dh_rxn and mass-balanced stoichiometry
        # assumption: this stoichiometry is mass-balanced
        self.stoich_init = Param(self.reaction_set, initialize={
                                            ("R1", "Vap", "H2O"): 4.95868,
                                            ("R1", "Vap", "CO2"): 5.95041556,
                                            ("R1", "Vap", "O2"): -6,
                                            ("R1", "Sol", "coal"): 0,
                                            ("R1", "Sol", "biomass"): -1,
                                            ("R1", "Vap", "N2"): 0,
                                            ("R1", "Sol", "ash"): 0.04409448,
                                                                    
                                            ("Rcoal", "Vap", "H2O"): 0.025,
                                            ("Rcoal", "Vap", "CO2"): 0.620971467,
                                            ("Rcoal", "Vap", "O2"): -0.548471467,
                                            ("Rcoal", "Sol", "coal"): -1,
                                            ("Rcoal", "Sol", "biomass"): 0,
                                            ("Rcoal", "Vap", "N2"): 0.005,
                                            ("Rcoal", "Sol", "ash"): 0.005361517,

                                            ("RCH4", "Vap", "H2O"):2,
                                            ("RCH4", "Vap", "CO2"):1,
                                            ("RCH4", "Vap", "O2"):-2,
                                            ("RCH4", "Vap", "CH4"):-1,
                                            }
                                            ,mutable=False)
        self.rate_reaction_stoichiometry = Var(self.reaction_set, initialize=self.stoich_init)
        self.rate_reaction_stoichiometry.fix()
        
        # self.reactant_list=Set(initialize=["biomass"])

        #fuel dict
        self.limit_reactant_dict = Param(self.rate_reaction_idx, initialize={
            "R1": ("Sol","biomass"),
            "Rcoal": ("Sol","coal"),
            "RCH4": ("Vap","CH4"),
        },
        within=Any)

        
            
            

        dh_rxn_dict = {"R1": -2749556.40, # @ w=9%, h=6% ==> ncv=-2749556.40
                       "Rcoal": -284675.1254, #[J/molCoal]
                       "RCH4": -802125 #based on 50 Mj/kg LHV
                       } 
        
        self.dh_rxn = Var(self.rate_reaction_idx, 
                          initialize = dh_rxn_dict,
                          domain=Reals,
                          doc="Heat of reaction")
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
    

