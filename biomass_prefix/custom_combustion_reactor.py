#custom_stoichiometric_reactor
#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Standard IDAES STOICHIOMETRIC reactor model
"""

# Import Pyomo libraries
from pyomo.environ import Reference, Var, Param, units as pyunits, value
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.models_extra.power_generation.unit_models.boiler_heat_exchanger import BoilerHeatExchanger
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block,
)

from biomass_combustion_rp import BMCombReactionParameterBlock



__author__ = "Chinedu Okoli, Andrew Lee"


@declare_process_block_class("MultiCombReactor")
class MultiCombReactorData(UnitModelBlockData):
    """
    Standard Stoichiometric Reactor Unit Model Class
    This model assumes that all given reactions are irreversible, and that each
    reaction has a fixed rate_reaction extent which has to be specified by the
    user.
    """

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.useDefault.
**Valid values:** {
**MaterialBalanceType.useDefault - refer to property package for default
balance type
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.useDefault,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.useDefault.
**Valid values:** {
**EnergyBalanceType.useDefault - refer to property package for default
balance type
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be constructed,
**default** - MomentumBalanceType.pressureTotal.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "has_heat_of_reaction",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Heat of reaction term construction flag",
            doc="""Indicates whether terms for heat of reaction terms should be
constructed,
**default** - False.
**Valid values:** {
**True** - include heat of reaction terms,
**False** - exclude heat of reaction terms.}""",
        ),
    )
    CONFIG.declare(
        "has_heat_transfer",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Heat transfer term construction flag",
            doc="""Indicates whether terms for heat transfer should be constructed,
**default** - False.
**Valid values:** {
**True** - include heat transfer terms,
**False** - exclude heat transfer terms.}""",
        ),
    )
    CONFIG.declare(
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )

#     CONFIG.declare(
#         "reaction_package",
#         ConfigValue(
#             default=None,
#             domain=is_reaction_parameter_block,
#             description="Reaction package to use for control volume",
#             doc="""Reaction parameter object used to define reaction calculations,
# **default** - None.
# **Valid values:** {
# **None** - no reaction package,
# **ReactionParameterBlock** - a ReactionParameterBlock object.}""",
#         ),
#     )
#     CONFIG.declare(
#         "reaction_package_args",
#         ConfigBlock(
#             implicit=True,
#             description="Arguments to use for constructing reaction packages",
#             doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
# and used when constructing these,
# **default** - None.
# **Valid values:** {
# see reaction package for documentation.}""",
#         ),
#     )


    def build(self):
        """
        Begin building model (pre-DAE transformation).
        Args:
            None
        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(MultiCombReactorData, self).build()
        self.reaction_package = BMCombReactionParameterBlock(property_package=self.config.property_package)
        # Build Control Volume
        self.control_volume = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
            reaction_package=self.reaction_package,
        )

        self.control_volume.add_state_blocks(has_phase_equilibrium=False)

        self.control_volume.add_reaction_blocks(has_equilibrium=False)

        self.control_volume.add_material_balances(
            balance_type=self.config.material_balance_type, has_rate_reactions=True
        )

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=self.config.has_heat_transfer,
            has_heat_of_reaction=self.config.has_heat_of_reaction,
        )

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )


        # Set references to balance terms at unit level
        if (
            self.config.has_heat_transfer is True
            and self.config.energy_balance_type != EnergyBalanceType.none
        ):
            self.heat_duty = Reference(self.control_volume.heat[:])
        if (
            self.config.has_pressure_change is True
            and self.config.momentum_balance_type != MomentumBalanceType.none
        ):
            self.deltaP = Reference(self.control_volume.deltaP[:])

        # Add Ports
        self.add_inlet_port()
        self.add_outlet_port()

        #creating conversion variable
        

        #add heat loss config option to turn off/on heat loss correlation. (may also need logic for front end state vars implementation?)
        #Q_loss = UAdT
        self.ohtc = Var(initialize=250, units=pyunits.J/pyunits.m**2/pyunits.K/pyunits.s)
        self.surface_area = Var(initialize=0.02, units=pyunits.m**2, doc="casing outer surface area")
        self.surface_temp = Var(initialize=55+273.12, units=pyunits.K, doc="outer skin temperature of boiler")
            
        self.hcon=Var(initialize=0.06) #concentration of hydrogen as a percentage of weight, h=6%
        self.wcon=Var(initialize=0.09) #water content of fuel as percentage of weight
        self.gcv=Param(initialize=20.2, units=pyunits.MJ/pyunits.kg, doc="gross calorific value") 

        


        #alternate build for conversion:
        # self.conversion = Var(self.reaction_package.rate_reaction_idx ,initialize=1, bounds=(0,1), units="dimensionless")
        for r in self.reaction_package.rate_reaction_idx:
            setattr(self,f"conversion_{r}", Var(initialize=1,bounds=(0,1), units="dimensionless"))
            setattr(self,f"dh_rxn_{r}", Var(initialize=1000000, units=pyunits.J/pyunits.mol))   
        

        for u in self.reaction_package.uncombs_set:
            setattr(self,f"ash_mass_{u}",Var(initialize=0.000001, units=pyunits.g/pyunits.g))

        #abbreviations
        mw = self.config.property_package.config.components
        # rrstoich = self.reaction_package.rate_reaction_stoichiometry #actual stoichiometry variables
        
        @self.Constraint(self.reaction_package.uncombs_set)
        def ash_con(b,u):
            p,l = b.reaction_package.limit_reactant_dict[u]
            ash_perc = getattr(b,f"ash_mass_{u}")
            ashi = b.reaction_package.stoich_init[u,"Sol","ash"]
            fueli = b.reaction_package.stoich_init[u,p,l]
            mw_ash = mw["ash"]["parameter_data"]["mw"][0]
            mw_fuel = mw[l]["parameter_data"]["mw"][0]
            b.reaction_package.rate_reaction_stoichiometry[u,"Sol","ash"].unfix()
            return b.reaction_package.rate_reaction_stoichiometry[u,"Sol","ash"]== (ash_perc*mw_fuel/mw_ash*(-fueli)/(1-ashi))*mw_ash/mw_fuel+ashi
        
        @self.Constraint(self.reaction_package.uncombs_set)
        def ash_con_fuel(b,u):
            p,l = self.reaction_package.limit_reactant_dict[u]
            ash_perc = getattr(b,f"ash_mass_{u}")
            ashi = b.reaction_package.stoich_init[u,"Sol","ash"] #initial ash is assumed to be part of the mass balance
            mw_ash = mw["ash"]["parameter_data"]["mw"][0]
            fueli = b.reaction_package.stoich_init[u,p,l]
            mw_fuel = mw[l]["parameter_data"]["mw"][0]
            b.reaction_package.rate_reaction_stoichiometry[u,p,l].unfix()
            return b.reaction_package.rate_reaction_stoichiometry[u,p,l] == -(ash_perc*mw_fuel/mw_ash*(-fueli)/(1-ashi))+fueli




        @self.Constraint(self.reaction_package.rate_reaction_idx)
        def dh_rxn_link(b,r):
            b.reaction_package.dh_rxn[r].unfix()
            return getattr(b,f"dh_rxn_{r}") == b.reaction_package.dh_rxn[r]

        #hard coded constraint just for biomassbiomass heating value [turn into callable method?]
        @self.Constraint()
        def ncv_eqn(b):
            p,l = self.reaction_package.limit_reactant_dict[u]
            ash_perc = getattr(b,f"ash_mass_{u}")
            mw_ash = mw["ash"]["parameter_data"]["mw"][0]
            ashi = b.reaction_package.stoich_init[u,"Sol","ash"]
            fueli = b.reaction_package.stoich_init[u,p,l]
            mw_fuel = mw[l]["parameter_data"]["mw"][0]
            return b.dh_rxn_R1 == (
                -(b.gcv*(1-b.wcon)-2.447*b.wcon-2.447*b.hcon*9.01*(1-b.wcon))*162.1394*1000/(-(ash_perc*mw_fuel/mw_ash*(-fueli)/(1-ashi))+fueli))
            
        
        @self.Constraint(self.flowsheet().time,)
        def heat_loss_eqn(b,t):
            return b.heat_duty[t] == (
            b.ohtc*b.surface_area*(-b.outlet.temperature[0]+b.surface_temp)
            )


        @self.Constraint(
                self.flowsheet().time,
                self.reaction_package.rate_reaction_idx)
        def conversion_performance_eqn(b, t, r):
            p,l = self.reaction_package.limit_reactant_dict[r]
            return getattr(b, f"conversion_{r}") == (
            b.control_volume.rate_reaction_extent[t,r]
             /(b.control_volume.properties_in[t].mole_frac_comp[l]
             *b.control_volume.properties_in[t].flow_mol
            ))
        


    def _get_performance_contents(self, time_point=0):
        var_dict = {
            "Water Content": self.wcon,
            "H2 Content": self.hcon,
            # "Biomass Calorific value": self.ncv,
            # "BM ncv": self.dh_rxn_Rbiomass
            # "overall heat transfer coefficient": self.ohtc,
            # "surface area": self.surface_area,
            # "surface temperature": self.surface_temp,
            }
        for r in self.reaction_package.rate_reaction_idx:
            var_dict["%s Conversion"%(r)] = getattr(self,f"conversion_{r}")
            var_dict["%s dh_rxn"%(r)] = getattr(self, f"dh_rxn_{r}")
        # for u in self.reaction_package.uncombs_set:
        #     var_dict["%s Ash content"%(u)] = getattr(self,f"ash_mass_{u}")
        if hasattr(self, "heat_duty"):
            var_dict["Heat Duty"] = self.heat_duty[time_point]
        if hasattr(self, "deltaP"):
            var_dict["Pressure Change"] = self.deltaP[time_point]

        return {"vars": var_dict}