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
        self.phase_list = Set(initialize=['Vap', 'Sol'])

        # Component list - a list of component identifiers
        self.component_list = Set(initialize=['H2O',
                                              'CO2',
                                              'O2',
                                              'CO',
                                              'N2',
                                              'biomass'])

        # Reaction Index
        self.rate_reaction_idx = Set(initialize=["R1"])

        self.ash_content = Var(initialize=0.01, doc="fraction of ash content in fuel. Assumes stoichiometry of fuel in reaction is -1 (see self.rate_reaction_stoichiometry)") 

        # Reaction Stoichiometry
        self.rate_reaction_stoichiometry = {("R1", "Vap", "H2O"): 5,
                                            ("R1", "Vap", "CO2"): 6,
                                            ("R1", "Vap", "O2"): -6,
                                            ("R1", "Sol", "biomass"): -1,
                                            ("R1", "Vap", "N2"): 0,
                                            ("R1", "Vap", "CO"): 0,
                                            ("R1", "Sol", "ash"): self.ash_content
                                            }

        
        self.reactant_list=Set(initialize=["biomass","O2"])
        
        self.h=Var(initialize=0.06) #concentration of hydrogen as a percentage of weight, h=6%
        self.w=Var(initialize=0.09) #water content of fuel as percentage of weight
        self.gcv=Param(initialize=20.2, units=pyunits.MJ/pyunits.kg, doc="gross calorific value") #gross calorific value (dry basis)
        self.ncv=(self.gcv*(1-self.w)-2.447*self.w-2.447*self.h*9.01*(1-self.w))*162.1394*1000 #J/mol 
        #net calorific value (wet basis) (pg. 7) https://www.mbie.govt.nz/dmsdocument/125-industrial-bioenergy-
        #ncv multiplied by 162 g/mo (cellulose) to convert from /mass to /mol basis.

        dh_rxn_dict = {"R1": -self.ncv} # w=9%, h=6% ncv=-2749556.40
        
        def dh_rxn(b,reaction_index):
            # only one reaction index, so we are just setting it to the ncv
            return -self.ncv

        self.dh_rxn = Expression(self.rate_reaction_idx, 
                            rule=dh_rxn,
                            doc="Heat of reaction")

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
    

