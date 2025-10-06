#custom_combustion_reactor

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

@declare_process_block_class("MultiCombReactor")
class MultiCombReactorData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,"""
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.useDefault,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,""",
        ),
    )
    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be constructed""",
        ),
    )
    CONFIG.declare(
        "has_heat_of_reaction",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Heat of reaction term construction flag",
            doc="""Indicates whether terms for heat of reaction terms should be constructed,""",
        ),
    )
    CONFIG.declare(
        "has_heat_transfer",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Heat transfer term construction flag",
            doc="""Indicates whether terms for heat transfer should be constructed,""",
        ),
    )
    CONFIG.declare(
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be constructed,""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)""",
        ),
    )

    def build(self):
        # Call UnitModel.build to setup dynamics
        super(MultiCombReactorData, self).build()
        self.reaction_package = BMCombReactionParameterBlock(property_package=self.config.property_package)

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

        #abbreviations
        mw = self.config.property_package.config.components        
        mw_air = 28.97
        O2_compound = "O2"
        N2_compound = "N2"

        # Add Custom Variables

        #casing heat loss variables
        self.ohtc = Var(initialize=250, units=pyunits.J/pyunits.m**2/pyunits.K/pyunits.s)
        self.surface_area = Var(initialize=0.02, units=pyunits.m**2, doc="casing outer surface area")
        self.surface_temp = Var(initialize=55+273.12, units=pyunits.K, doc="outer skin temperature of boiler")
            
        # biomass-specific NCV variables
        self.hcon=Var(initialize=0.06) #concentration of hydrogen as a percentage of weight, h=6%
        self.wcon=Var(initialize=0.09) #water content of fuel as percentage of weight
        self.gcv=Param(initialize=20.2, units=pyunits.MJ/pyunits.kg, doc="gross calorific value") 

        # variables for each reaction in rate_reaction_idx list
        for r in self.reaction_package.rate_reaction_idx:
            p,l = self.reaction_package.limit_reactant_dict[r]
            setattr(self,f"conversion_{r}", Var(initialize=1,bounds=(0,1), units="dimensionless"))
            setattr(self,f"dh_rxn_{r}", Var(initialize=1000000, units=pyunits.J/pyunits.mol))   
            setattr(self,f"excess_air_percent_{r}", Var(initialize = 15, units="dimensionless"))
            setattr(self,f"air_fuel_mass_ratio_{r}", Var(initialize=6, units="dimensionless"))
            setattr(self,f"mass_flow_kg_{l}", Var(initialize=1 ,doc="mass flow of fuel in kg/s"))
            setattr(self,f"mole_flow_{l}", Var(initialize=100, doc="fuel molar flow in mol/s"))
            setattr(self,f"mole_flow_air_{r}", Var(initialize=1000, doc="mole flow air in mol/s"))

        # ash mass content variable only for reactions with ash in uncombs_set list        
        for u in self.reaction_package.uncombs_set:
            setattr(self,f"ash_mass_{u}",Var(initialize=0, units=pyunits.g/pyunits.g))

        # Add Custom Constraints

        @self.Constraint(self.reaction_package.rate_reaction_idx)
        def fuel_flow_conversion(b,r): #kg/s to mol/s
            p,l = b.reaction_package.limit_reactant_dict[r]
            return getattr(b, f"mass_flow_kg_{l}") == getattr(b, f"mole_flow_{l}")*mw[l]["parameter_data"]["mw"][0]/1000
        
        @self.Constraint(self.flowsheet().time,self.reaction_package.rate_reaction_idx)
        def mole_flow_link(b,t,r):
            p,l = b.reaction_package.limit_reactant_dict[r]
            return getattr(b, f"mole_flow_{l}") == b.control_volume.properties_in[t].mole_frac_comp[l]*b.control_volume.properties_in[t].flow_mol
        
        @self.Constraint(self.reaction_package.rate_reaction_idx)
        def Air_Fuel_ratio(b,r):
            p,l = b.reaction_package.limit_reactant_dict[r]
            return getattr(b, f"air_fuel_mass_ratio_{r}") == getattr(b, f"mole_flow_air_{r}")*mw_air/(getattr(b, f"mass_flow_kg_{l}")*1000)
        
        @self.Constraint(self.reaction_package.rate_reaction_idx)
        def excess_air(b,r):
            p,l = b.reaction_package.limit_reactant_dict[r]
            O2_stoich = b.reaction_package.rate_reaction_stoichiometry[r,"Vap",O2_compound]
            fuel_stoich = b.reaction_package.rate_reaction_stoichiometry[r,p,l]
            return (1+getattr(b, f"excess_air_percent_{r}")/100)*(-O2_stoich)/0.21/(-fuel_stoich)*getattr(b, f"mole_flow_{l}") == getattr(b, f"mole_flow_air_{r}")
        
        @self.Constraint(self.flowsheet().time)
        def N2_flow_link(b,t):
            return sum(getattr(b, f"mole_flow_air_{r}") for r in b.reaction_package.rate_reaction_idx)*0.79 == b.control_volume.properties_in[t].mole_frac_comp[N2_compound]*b.control_volume.properties_in[t].flow_mol

        @self.Constraint(self.flowsheet().time)
        def O2_flow_link(b,t):
            return sum(getattr(b, f"mole_flow_air_{r}") for r in b.reaction_package.rate_reaction_idx)*0.21 == b.control_volume.properties_in[t].mole_frac_comp[O2_compound]*b.control_volume.properties_in[t].flow_mol

        @self.Constraint(self.flowsheet().time)
        def total_flow_link(b,t):
            l = b.reaction_package.limit_reactant_dict
            rxns = b.reaction_package.rate_reaction_idx
            return b.control_volume.properties_in[t].flow_mol == sum(getattr(b, f"mole_flow_{l[r][1]}") for r in rxns)+sum(getattr(b, f"mole_flow_air_{r}") for r in rxns)


        @self.Constraint(self.reaction_package.uncombs_set)
        def ash_con(b,u):
            p,l = b.reaction_package.limit_reactant_dict[u]
            ash_perc = getattr(b,f"ash_mass_{u}")
            ashi = b.reaction_package.stoich_init[u,"Sol","ash"]
            fueli = b.reaction_package.stoich_init[u,p,l]
            mw_ash = mw["ash"]["parameter_data"]["mw"][0]
            mw_fuel = mw[l]["parameter_data"]["mw"][0]
            ash_perc_mol = ash_perc*mw_fuel/mw_ash
            added_mols_BM = (ash_perc_mol-ashi)*(-fueli)/(1-(ash_perc_mol-ashi))
            b.reaction_package.rate_reaction_stoichiometry[u,"Sol","ash"].unfix() 
            return b.reaction_package.rate_reaction_stoichiometry[u,"Sol","ash"]== added_mols_BM*(1+ashi)+ashi*(-fueli)
        
        @self.Constraint(self.reaction_package.uncombs_set)
        def ash_con_fuel(b,u):
            p,l = self.reaction_package.limit_reactant_dict[u]
            ash_perc = getattr(b,f"ash_mass_{u}")
            ashi = b.reaction_package.stoich_init[u,"Sol","ash"] #initial ash is assumed to be part of the mass balance
            mw_ash = mw["ash"]["parameter_data"]["mw"][0]
            fueli = b.reaction_package.stoich_init[u,p,l]
            mw_fuel = mw[l]["parameter_data"]["mw"][0]
            ash_perc_mol = ash_perc*mw_fuel/mw_ash
            added_mols_BM = (ash_perc_mol-ashi)*(-fueli)/(1-(ash_perc_mol-ashi))
            b.reaction_package.rate_reaction_stoichiometry[u,p,l].unfix()
            return b.reaction_package.rate_reaction_stoichiometry[u,p,l] == -(added_mols_BM)+fueli

        @self.Constraint(self.reaction_package.rate_reaction_idx)
        def dh_rxn_link(b,r):
            b.reaction_package.dh_rxn[r].unfix()
            return getattr(b,f"dh_rxn_{r}") == b.reaction_package.dh_rxn[r]

        #dedicated NCV constaint for Biomass
        @self.Constraint()
        def ncv_eqn(b):
            ash_perc = getattr(b,f"ash_mass_R1")
            mw_ash = mw["ash"]["parameter_data"]["mw"][0]
            ashi = b.reaction_package.stoich_init["R1","Sol","ash"]
            fueli = b.reaction_package.stoich_init["R1","Sol","biomass"]
            mw_fuel = mw["biomass"]["parameter_data"]["mw"][0]
            ash_perc_mol = ash_perc*mw_fuel/mw_ash
            added_mols_BM = (ash_perc_mol-ashi)*(-fueli)/(1-(ash_perc_mol-ashi))
            return b.dh_rxn_R1 == (
                -(b.gcv*(1-b.wcon)-2.447*b.wcon-2.447*b.hcon*9.01*(1-b.wcon))*162.1394*1000*(-fueli)/(added_mols_BM-fueli))
            
        @self.Constraint(self.flowsheet().time,)
        def casing_heat_loss(b,t):
            return b.heat_duty[t] == (
            b.ohtc*b.surface_area*(-b.outlet.temperature[0]+b.surface_temp)
            )

        @self.Constraint(self.flowsheet().time,self.reaction_package.rate_reaction_idx)
        def conversion_performance_eqn(b, t, r):
            p,l = self.reaction_package.limit_reactant_dict[r]
            return getattr(b, f"conversion_{r}") == (
            b.control_volume.rate_reaction_extent[t,r]
             /(b.control_volume.properties_in[t].mole_frac_comp[l]
             *b.control_volume.properties_in[t].flow_mol
            ))

    #variables displayed in terminal unit model report
    def _get_performance_contents(self, time_point=0):
        var_dict = {
            "Biomass Water Content": self.wcon,
            "Biomass H2 Content": self.hcon,
            }
        for r in self.reaction_package.rate_reaction_idx:
            var_dict["Conversion %s"%(r)] = getattr(self,f"conversion_{r}")
            var_dict["dh_rxn %s"%(r)] = getattr(self, f"dh_rxn_{r}")
        for u in self.reaction_package.uncombs_set:
            var_dict["Ash content %s"%(u)] = getattr(self,f"ash_mass_{u}")
        if hasattr(self, "heat_duty"):
            var_dict["Heat Duty"] = self.heat_duty[time_point]
        if hasattr(self, "deltaP"):
            var_dict["Pressure Change"] = self.deltaP[time_point]
        return {"vars": var_dict}