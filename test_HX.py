# test HX

# Import Pyomo libraries
from pyomo.environ import ConcreteModel, SolverFactory, value
# Import IDAES core
from idaes.core import FlowsheetBlock
# Import Unit Model Modules
from idaes.models.properties import iapws95
# import ideal flue gas prop pack
# from test_pp import configuration
from idaes.models_extra.power_generation.properties.flue_gas_ideal import FlueGasParameterBlock
from adapted_flue_gas_ideal import FlueGasParameterBlock
from idaes.models_extra.gas_solid_contactors.properties.methane_iron_OC_reduction.solid_phase_thermo import SolidPhaseParameterBlock

from test_pp import configuration
from idaes.models.properties.modular_properties import GenericParameterBlock

# Import Power Plant HX Unit Model
from idaes.models_extra.power_generation.unit_models.boiler_heat_exchanger import (
    BoilerHeatExchanger,
    TubeArrangement,
    HeatExchangerFlowPattern,
)
from idaes.models.unit_models import (
Mixer,
StoichiometricReactor,
Heater,
HeatExchanger1D,
HeatExchanger,
Separator,
Flash
)
import pyomo.environ as pe # Pyomo environment
from pyomo.environ import units as pyunits
from idaes.core import FlowsheetBlock, StateBlock
from idaes.models.unit_models.heat_exchanger import delta_temperature_amtd_callback
from idaes.models.properties.modular_properties.pure.ChapmanEnskog import collision_integral_neufeld_callback 

from idaes.models.properties import iapws95

from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.core.util.model_diagnostics import (
    DiagnosticsToolbox,
)
# Create a Concrete Model as the top level object
m = ConcreteModel()

# Add a flowsheet object to the model
m.fs = FlowsheetBlock(dynamic=False)

# Add property packages to flowsheet library
m.fs.steam_properties = iapws95.Iapws95ParameterBlock()
m.fs.biomass_properties = GenericParameterBlock(**configuration)
# m.fs.biomass_properties = FlueGasParameterBlock(components=["N2", "O2", "CO2", "H2O","ash", "biomass"])

# Create unit models
m.fs.ECON = BoilerHeatExchanger(
    cold_side_name="tube",
    hot_side_name="shell",
    shell={"property_package": m.fs.biomass_properties,
           "has_pressure_change": True,
           },
    tube={"property_package": m.fs.steam_properties,
          "has_pressure_change": True,
          },
    delta_temperature_callback=delta_temperature_amtd_callback,
    # has_radiation=False
)


h = iapws95.htpx(463.706*pyunits.K, 2.5449e7*pyunits.Pa)
h2 = iapws95.htpx(523.706*pyunits.K, 2.5449e7*pyunits.Pa)

m.fs.ECON.tube_inlet.flow_mol[0].fix(24678.26) # mol/s
m.fs.ECON.tube_inlet.enth_mol[0].fix(h)
m.fs.ECON.tube_inlet.pressure[0].fix(2.5449e7) # Pa
m.fs.ECON.tube_outlet.enth_mol[0].fix(h2)

# FLUE GAS Inlet from Primary Superheater
FGrate = 60.3876e2  # mol/s equivalent of ~1930.08 klb/hr

m.fs.ECON.shell_inlet.mole_frac_comp[0,"H2O"].fix(0.4)
m.fs.ECON.shell_inlet.mole_frac_comp[0,"O2"].fix(0.4)
m.fs.ECON.shell_inlet.mole_frac_comp[0,"N2"].fix(0)
m.fs.ECON.shell_inlet.mole_frac_comp[0,"CO2"].fix(0)
m.fs.ECON.shell_inlet.mole_frac_comp[0,"ash"].fix(0.1)
m.fs.ECON.shell_inlet.mole_frac_comp[0,"biomass"].fix(0.1)
m.fs.ECON.shell_inlet.temperature[0].fix(682.335) # K
m.fs.ECON.shell_inlet.pressure[0].fix(100145) # Pa
m.fs.ECON.shell_inlet.flow_mol[0].fix(FGrate) # Pa


# m.fs.ECON.area.fix(30000)
# m.fs.ECON.overall_heat_transfer_coefficient.fix(170)
# m.fs.ECON.tube_outlet.pressure[0].fix(2.4449e7)
# m.fs.ECON.shell_outlet.pressure[0].fix(99145)

#1D specific:
ITM = 0.0254  # inch to meter conversion
m.fs.ECON.tube_di.fix((2-2*0.188)*ITM)
m.fs.ECON.tube_thickness.fix(0.188*ITM)
m.fs.ECON.pitch_x.fix(3.5*ITM)
# pitch_y = (54.5) gas path transverse width /columns
m.fs.ECON.pitch_y.fix(5.03*ITM)
m.fs.ECON.tube_length.fix(53.41*12*ITM) # use tube length (53.41 ft)
m.fs.ECON.tube_nrow.fix(36*2.5)         # use to match baseline performance
m.fs.ECON.tube_ncol.fix(130)            # 130 from NETL report
m.fs.ECON.nrow_inlet.fix(2)
m.fs.ECON.delta_elevation.fix(50)

m.fs.ECON.tube_rfouling = 0.000176
m.fs.ECON.shell_rfouling = 0.00088
if m.fs.ECON.config.has_radiation is True:
    m.fs.ECON.emissivity_wall.fix(0.7)       

m.fs.ECON.fcorrection_htc.fix(1.5)
m.fs.ECON.fcorrection_dp_tube.fix(1.0)
m.fs.ECON.fcorrection_dp_shell.fix(1.0)

# Initialize the model
print(degrees_of_freedom(m))  
m.fs.ECON.initialize()

m.fs.ECON.report()

dt = DiagnosticsToolbox(m)
dt.report_structural_issues()
dt.display_underconstrained_set()