from pyomo.environ import (
    Constraint,             # 约束
    Var,                    # 变量
    ConcreteModel,          # 模型
    Expression,             # 表达式
    Objective,              # 目标
    SolverFactory,          # 求解器
    TransformationFactory,  # 转换器，应用模型变换（如离散化）
    value,                  # 取出值
)
from idaes.core import FlowsheetBlock                    # 流程图
from pyomo.network import Arc, SequentialDecomposition   # 连接单元模型； 处理循环

import props as props
import props_reaction as props_reaction
from idaes.models.properties.iapws95 import Iapws95ParameterBlock

from idaes.models.properties.modular_properties.base.generic_property import \
    GenericParameterBlock
from idaes.models.properties.modular_properties.base.generic_reaction import \
    GenericReactionParameterBlock
from idaes.models.unit_models import (
    Feed,                               # 进料
    PressureChanger,                    # 压力机
    Mixer,                              # 混合器
    Separator,                          # 分离器
    Separator as Splitter,              # 分离器
    SplittingType,                      # 分离器
    Heater,                             # 加热器
    StoichiometricReactor,              # 反应器
    Flash,                              # （闪蒸）蒸馏器
    HeatExchanger,                      # 热交换器
    Feed,                               # 进料口
    EquilibriumReactor,                 # 平衡反应器
    Product                             # 出口
)
from idaes.models.unit_models.heat_exchanger import (
    HeatExchanger,
    delta_temperature_amtd_callback,
)
from idaes.models.unit_models.pressure_changer import (
    ThermodynamicAssumption,

)
from idaes.core.util.model_statistics import degrees_of_freedom

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

m.fs.props_water = Iapws95ParameterBlock()

m.fs.props = GenericParameterBlock(**props.configuration)

m.fs.props_reaction = GenericReactionParameterBlock(
    property_package=m.fs.props,
    **props_reaction.config_dict)



# 电解槽后分理出氢气流氧气流
m.fs.Electrolyzer = StoichiometricReactor(
    property_package=m.fs.props,
    reaction_package=m.fs.props_reaction,
    has_heat_of_reaction=True,
    has_heat_transfer=True,
    has_pressure_change=False,
)
m.fs.Sep = Separator(
    property_package=m.fs.props,
    outlet_list=["H2_outlet", "O2_outlet"],
    split_basis=SplittingType.componentFlow
)


# 氢分离器，要定义夹带率
m.fs.Flash_H2_1 = Flash(
    property_package=m.fs.props,
    has_heat_transfer=True,
    has_pressure_change=True
)
m.fs.Split_liq_1 = Splitter(
    property_package=m.fs.props,
    outlet_list=["to_recovery", "to_mixer"]
)
m.fs.Mixer_H2_1 = Mixer(
    property_package=m.fs.props,
    inlet_list=["liq_1", "flash_vap"],
    )




# 氢气换热器
m.fs.HE_H2 = HeatExchanger(
    delta_temperature_callback=delta_temperature_amtd_callback,
    hot_side_name="shell",
    cold_side_name="tube",
    shell={"property_package": m.fs.props},
    tube={"property_package": m.fs.props_water},
)


# 氢洗涤器
m.fs.Mixer_H2_2_1 = Mixer(
    property_package=m.fs.props,
    inlet_list=["KOH_Feed", "H2_from_HE"],
)
m.fs.Flash_H2_2 = Flash(
    property_package=m.fs.props,
    has_heat_transfer=True,
    has_pressure_change=True
    )
m.fs.Split_liq_2 = Splitter(
    property_package=m.fs.props,
    outlet_list=["to_recovery", "to_mixer"]
)
m.fs.Mixer_H2_2_2 = Mixer(
    property_package=m.fs.props,
    inlet_list=["liq_2", "flash_vap"],
    )




# 供冷却水和除盐水
m.fs.H2O1 = Feed(property_package=m.fs.props_water)
m.fs.H2O2 = Feed(property_package=m.fs.props_water)
m.fs.KOH = Feed(property_package=m.fs.props)

# 碱液回收、换热、加压
m.fs.Mixer_H2_3 = Mixer(
    property_package=m.fs.props,
    inlet_list=["KOH_outlet1", "KOH_outlet2"],
    )
m.fs.HE_KOH = HeatExchanger(
    delta_temperature_callback=delta_temperature_amtd_callback,
    hot_side_name="shell",
    cold_side_name="tube",
    shell={"property_package": m.fs.props},
    tube={"property_package": m.fs.props_water},
    )
m.fs.Press_KOH = PressureChanger(
    property_package=m.fs.props,
    thermodynamic_assumption=ThermodynamicAssumption.isothermal,
    has_pressure_change=True
    )



# 暂定为出口
m.fs.H2_product = Product(property_package=m.fs.props)
m.fs.O2_product = Product(property_package=m.fs.props)
m.fs.CW_out1     = Product(property_package=m.fs.props_water)
m.fs.CW_out2     = Product(property_package=m.fs.props_water)




m.fs.s01 = Arc(source=m.fs.Electrolyzer.outlet, destination=m.fs.Sep.inlet)
m.fs.s02 = Arc(source=m.fs.Sep.H2_outlet, destination=m.fs.Flash_H2_1.inlet)
m.fs.s03 = Arc(source=m.fs.Sep.O2_outlet, destination=m.fs.O2_product)

m.fs.s04 = Arc(source=m.fs.Flash_H2_1.liq_outlet, destination=m.fs.Split_liq_1.inlet)
m.fs.s05 = Arc(source=m.fs.Split_liq_1.to_mixer, destination=m.fs.Mixer_H2_1.liq_1)
m.fs.s06 = Arc(source=m.fs.Flash_H2_1.vap_outlet, destination=m.fs.Mixer_H2_1.flash_vap)

m.fs.s07 = Arc(source=m.fs.Mixer_H2_1.outlet, destination=m.fs.HE_H2.shell_inlet)
m.fs.s08 = Arc(source=m.fs.H2O1.outlet, destination=m.fs.HE_H2.tube_inlet)
m.fs.s023 = Arc(source=m.fs.HE_H2.tube_outlet, destination=m.fs.CW_out1)

m.fs.s09 = Arc(source=m.fs.HE_H2.shell_outlet, destination=m.fs.Mixer_H2_2_1.H2_from_HE)
m.fs.s10 = Arc(source=m.fs.KOH.outlet, destination=m.fs.Mixer_H2_2_1.KOH_Feed)

m.fs.s11 = Arc(source=m.fs.Mixer_H2_2_1.outlet, destination=m.fs.Flash_H2_2.inlet)
m.fs.s12 = Arc(source=m.fs.Flash_H2_2.liq_outlet, destination=m.fs.Split_liq_2.inlet)
m.fs.s13 = Arc(source=m.fs.Split_liq_2.to_mixer, destination=m.fs.Mixer_H2_2_2.liq_2)
m.fs.s14 = Arc(source=m.fs.Flash_H2_2.vap_outlet, destination=m.fs.Mixer_H2_2_2.flash_vap)

m.fs.s15 = Arc(source=m.fs.Mixer_H2_2_2.outlet, destination=m.fs.H2_product)

m.fs.s16 = Arc(source=m.fs.Split_liq_1.to_recovery, destination=m.fs.Mixer_H2_3.KOH_outlet1)
m.fs.s17 = Arc(source=m.fs.Split_liq_2.to_recovery, destination=m.fs.Mixer_H2_3.KOH_outlet2)

m.fs.s18 = Arc(source=m.fs.Mixer_H2_3.outlet, destination=m.fs.HE_KOH.shell_inlet)
m.fs.s19 = Arc(source=m.fs.H2O2.outlet, destination=m.fs.HE_KOH.tube_inlet)
m.fs.s024 = Arc(source=m.fs.HE_KOH.tube_outlet, destination=m.fs.CW_out2)

m.fs.s20 = Arc(source=m.fs.HE_KOH.shell_outlet, destination=m.fs.Press_KOH.inlet)
m.fs.s21 = Arc(source=m.fs.Press_KOH.outlet, destination=m.fs.Electrolyzer.inlet)



TransformationFactory("network.expand_arcs").apply_to(m)    # 将连接转换为变量
print()
print(degrees_of_freedom(m))





