from pyomo.environ import (
    ConcreteModel,
    Var,
    Constraint,
    TransformationFactory,
    value,
    Block,
    Expression
)

from idaes.core import FlowsheetBlock
from pyomo.network import Arc
from copy import deepcopy

import idaes.logger as idaeslog
from idaes.core.solvers import get_solver

import props as props
import props_reaction as props_reaction
import rxn_cat as rxn_cat
import rxn_an as rxn_an

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
)
from idaes.models.properties.modular_properties.state_definitions import FpcTP
from idaes.models.unit_models import (
    Feed,
    Mixer,
    Separator,
    SplittingType,
    StoichiometricReactor,
    Flash,
    Product,
    Heater,
    PressureChanger,
    Separator as Splitter,  # 分离器
    Translator,
)

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
import time, os
from idaes.core.util.exceptions import InitializationError     # 初始化错误
from idaes.core.util.model_serializer import to_json, from_json
from Electronchemical0D import AlkalineElectrolyzer0D
from creating_model import GasWasher
from translator import CompositionTranslator
import Liq as liq
import Vap as gas







# Build model
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

# Property and reaction packages
m.fs.props = GenericParameterBlock(**props.configuration)
m.fs.gas_props = GenericParameterBlock(**gas.configuration)
m.fs.liq_props = GenericParameterBlock(**liq.configuration)
m.fs.rxn = GenericReactionParameterBlock(property_package=m.fs.props, **props_reaction.config_dict)

m.fs.ely = AlkalineElectrolyzer0D(
    property_package=m.fs.props,
    reaction_package=m.fs.rxn
)
m.fs.Sep_0 = Separator(
    property_package=m.fs.props,
    outlet_list=["h2", "o2"],
    split_basis=SplittingType.phaseComponentFlow,
    )
m.fs.Sep_1 = Separator(
    property_package=m.fs.props,
    outlet_list=["to_washer", "to_splitter"],
    split_basis=SplittingType.phaseComponentFlow,
    )
m.fs.Sep_2 = Separator(
    property_package=m.fs.props,
    outlet_list=["to_recycle", "to_mix"],
    split_basis=SplittingType.phaseComponentFlow,
    )

m.fs.gin = CompositionTranslator(
    inlet_property_package=m.fs.props,
    outlet_property_package=m.fs.gas_props,
    mapping={
        ("Vap", "H2"): ("Vap", "H2"),
        ("Vap", "O2"): ("Vap", "O2"),
        ("Vap", "H2O"): ("Vap", "H2O"),
    }
)
m.fs.lin = CompositionTranslator(
    inlet_property_package=m.fs.props,
    outlet_property_package=m.fs.liq_props,
    mapping={
        ("Liq","H2O"): ("Liq","H2O"),
        ("Liq","KOH"): ("Liq","KOH"),
    }
)
m.fs.gout = CompositionTranslator(
    inlet_property_package=m.fs.gas_props,
    outlet_property_package=m.fs.props,
    mapping={
        ("Vap", "H2"): ("Vap", "H2"),
        ("Vap", "O2"): ("Vap", "O2"),
        ("Vap", "H2O"): ("Vap", "H2O"),
    }
)
m.fs.lout = CompositionTranslator(
    inlet_property_package=m.fs.liq_props,
    outlet_property_package=m.fs.props,
    mapping={
        ("Liq","H2O"): ("Liq","H2O"),
        ("Liq","KOH"): ("Liq","KOH"),
    }
)
m.fs.washer = GasWasher(
    gas_property_package=m.fs.gas_props,
    liq_property_package=m.fs.liq_props,
)
m.fs.mixer = Mixer(
    property_package=m.fs.props, momentum_mixing_type=MomentumMixingType.equality,
    inlet_list=["liq", "gas"],
    )

m.fs.water = Feed(property_package=m.fs.props)

m.fs.O2_product = Product(property_package=m.fs.props)
m.fs.H2_product = Product(property_package=m.fs.props)
m.fs.recycle = Product(property_package=m.fs.props)


m.fs.a01 = Arc(source=m.fs.ely.outlet, destination=m.fs.Sep_0.inlet)

m.fs.a02 = Arc(source=m.fs.Sep_0.o2, destination=m.fs.O2_product.inlet)
m.fs.a03 = Arc(source=m.fs.Sep_0.h2, destination=m.fs.Sep_1.inlet)

m.fs.a04 = Arc(source=m.fs.Sep_1.to_washer, destination=m.fs.gin.inlet)
m.fs.a05 = Arc(source=m.fs.gin.outlet, destination=m.fs.washer.gas_inlet)
m.fs.a06 = Arc(source=m.fs.Sep_1.to_splitter, destination=m.fs.Sep_2.inlet)

m.fs.a07 = Arc(source=m.fs.Sep_2.to_recycle, destination=m.fs.recycle.inlet)
m.fs.a08 = Arc(source=m.fs.Sep_2.to_mix, destination=m.fs.mixer.gas)

m.fs.a09 = Arc(source=m.fs.water.outlet, destination=m.fs.lin.inlet)
m.fs.a10 = Arc(source=m.fs.lin.outlet, destination=m.fs.washer.liquid_inlet)
m.fs.a11 = Arc(source=m.fs.washer.gas_outlet, destination=m.fs.gout.inlet)
m.fs.a12 = Arc(source=m.fs.gout.outlet, destination=m.fs.mixer.liq)

m.fs.a13 = Arc(source=m.fs.washer.liquid_outlet, destination=m.fs.lout.inlet)
m.fs.a14 = Arc(source=m.fs.lout.outlet, destination=m.fs.recycle.inlet)

m.fs.a15 = Arc(source=m.fs.mixer.outlet, destination=m.fs.H2_product.inlet)

TransformationFactory("network.expand_arcs").apply_to(m.fs)

t = m.fs.time.first()
print("Initial DOF =", degrees_of_freedom(m))

# 电解槽参数
m.fs.ely.inlet.temperature.fix(343.15)
m.fs.ely.inlet.pressure.fix(16e5)
m.fs.ely.current.fix(5000)
m.fs.ely.volume.fix(1.0)

# 电解槽进料
m.fs.ely.inlet.flow_mol_phase_comp[t, "Liq", "H2O"].fix(9)
m.fs.ely.inlet.flow_mol_phase_comp[t, "Liq", "KOH"].fix(1)
m.fs.ely.inlet.flow_mol_phase_comp[t, "Vap", "H2"].fix(1e-6)
m.fs.ely.inlet.flow_mol_phase_comp[t, "Vap", "O2"].fix(1e-6)
m.fs.ely.inlet.flow_mol_phase_comp[t, "Vap", "H2O"].fix(1e-6)







# m.fs.Flash_A = Flash(
#     property_package=m.fs.props, has_heat_transfer=True, has_pressure_change=True
# )


# # m.fs.mixer.feed.flow_mol_phase_comp[0, "Vap", "H2"].fix(1e-5)
# # m.fs.mixer.feed.flow_mol_phase_comp[0, "Vap", "O2"].fix(1e-5)
# # m.fs.mixer.feed.flow_mol_phase_comp[0, "Vap", "H2O"].fix(1e-5)
# # m.fs.mixer.feed.flow_mol_phase_comp[0, "Liq", "H2O"].fix(1)
# # m.fs.mixer.feed.flow_mol_phase_comp[0, "Liq","OH_ion"].fix(0.1)
# # m.fs.mixer.feed.flow_mol_phase_comp[0, "Liq", "K_ion"].fix(0.1)
# # m.fs.mixer.feed.temperature.fix(303.2)
# # m.fs.mixer.feed.pressure.fix(350000)
# m.fs.Anode.inlet.flow_mol_phase_comp[0, "Vap", "H2"].fix(1e-5)
# m.fs.Anode.inlet.flow_mol_phase_comp[0, "Vap", "O2"].fix(1e-5)
# m.fs.Anode.inlet.flow_mol_phase_comp[0, "Vap", "H2O"].fix(1e-5)
# m.fs.Anode.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(8)
# m.fs.Anode.inlet.flow_mol_phase_comp[0, "Liq","KOH"].fix(1)
# m.fs.Anode.inlet.temperature.fix(303.2)
# m.fs.Anode.inlet.pressure.fix(1900000)
#
#
#
# # 做一个简化，把I变成反应的摩尔流量吧
# # m.fs.F_const = 96485.0# C/mol
# m.fs.I = Var(initialize=0.5)
# m.fs.I.fix(0.5)
# m.fs.F_const = 0.25
# m.fs.anode_faraday = Constraint(
#     expr=m.fs.Anode.control_volume.rate_reaction_extent[0.0, "anode"]
#          == m.fs.I / (4 * m.fs.F_const)
# )
# m.fs.Anode.heat_duty.fix(0)
# # m.fs.Anode.outlet.temperature[0].fix(343.2)
#
#
# m.fs.eq_flash_temp = Constraint(
#     expr=m.fs.Flash_A.control_volume.properties_out[0].temperature
#          == m.fs.Anode.outlet.temperature[0]
# )
# m.fs.Flash_A.deltaP.fix(-1e6)
#
# for j in ["H2", "O2", "H2O"]:
#     m.fs.Spi.split_fraction[0, "to_mix","Vap", j].fix(1e-7)
# for j in ["K_ion", "H2O","OH_ion"]:
#     m.fs.Spi.split_fraction[0, "to_mix","Liq",  j].fix(0.1)
# # m.fs.Sep_A1.split_fraction[0, "to_O2", "Vap", "O2"].fix(1.0)
# # m.fs.Sep_A1.split_fraction[0, "to_O2", "Vap", "H2"].fix(1e-4)
# # m.fs.Sep_A1.split_fraction[0, "to_O2", "Vap", "H2O"].fix(1e-4)
# # m.fs.Sep_A1.split_fraction[0, "to_O2", "Liq", "H2O"].fix(1e-4)
# # m.fs.Sep_A1.split_fraction[0, "to_O2", "Liq", "OH_ion"].fix(1e-4)
# # m.fs.Sep_A1.split_fraction[0, "to_O2", "Liq", "K_ion"].fix(1e-4)
# print("Initial DOF =", degrees_of_freedom(m))
# # print("Anode DOF =", degrees_of_freedom(m.fs.Spi))
#
#
#
# m.fs.Anode.initialize(outlvl=idaeslog.INFO_LOW)
# # m.fs.Flash_A.initialize(outlvl=idaeslog.INFO_LOW)
# # m.fs.Sep_A.initialize(outlvl=idaeslog.INFO_LOW)
# # m.fs.Sep_A1.initialize(outlvl=idaeslog.INFO_LOW)
# # m.fs.mixer.initialize(outlvl=idaeslog.INFO_LOW)
# m.fs.O2_product.initialize(outlvl=idaeslog.INFO_LOW)
#
#
# solver = get_solver()
# solver.options["max_iter"] = 5000
# results = solver.solve(m, tee=True)

if os.path.exists("Hydrodealkylation.json"):
    os.remove("Hydrodealkylation.json")
m.fs.visualize("Hydrodealkylation")
while True:
    time.sleep(1)