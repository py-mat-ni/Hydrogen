from pyomo.environ import (
    ConcreteModel,
    Var,
    Constraint,
    ConstraintList,
    TransformationFactory,
    value,
    Block,
    Expression
)
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor
from idaes.core.util.scaling import constraint_scaling_transform
from idaes.core.util import DiagnosticsToolbox
import pyomo.environ as pyo
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
import pandas as pd
def add_scaling(m):
    # # 关键变量缩放
    # if hasattr(m.fs.ely, "current"): set_scaling_factor(m.fs.ely.current, 1e-4)
    # if hasattr(m.fs.ely, "voltage"): set_scaling_factor(m.fs.ely.voltage, 1e0)
    # if hasattr(m.fs.ely, "power"): set_scaling_factor(m.fs.ely.power, 1e-6)
    # if hasattr(m.fs.ely, "Pv_KOH"): set_scaling_factor(m.fs.ely.Pv_KOH, 1e-5)

    # 全局流量与焓值
    for v in m.component_data_objects(Var, descend_into=True):
        if "flow_mol_phase_comp" in v.name: set_scaling_factor(v, 1e-2)
        # if "enth_mol_phase" in v.name: set_scaling_factor(v, 1e-4)

    calculate_scaling_factors(m)
def build_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # 1. 导入物性包
    m.fs.props = GenericParameterBlock(**props.configuration_VLE)
    m.fs.props0 = GenericParameterBlock(**props.configuration_ELY)
    m.fs.rxn = GenericReactionParameterBlock(property_package=m.fs.props0, **props_reaction.config_dict)

    # 2. 单元定义
    m.fs.ely = AlkalineElectrolyzer0D(property_package=m.fs.props0, reaction_package=m.fs.rxn)
    m.fs.translator = Translator(inlet_property_package=m.fs.props0, outlet_property_package=m.fs.props)
    m.fs.translator1 = Translator(inlet_property_package=m.fs.props, outlet_property_package=m.fs.props0)

    m.fs.Sep_0 = Separator(property_package=m.fs.props, outlet_list=["h2", "o2"],
                           split_basis=SplittingType.phaseComponentFlow)
    m.fs.Sep_O2 = Separator(property_package=m.fs.props, outlet_list=["recycle", "o2"],
                            split_basis=SplittingType.phaseComponentFlow)
    m.fs.Sep_H2 = Separator(property_package=m.fs.props, outlet_list=["h2", "recycle"],
                            split_basis=SplittingType.phaseComponentFlow)

    m.fs.Heater_O2 = Heater(property_package=m.fs.props, has_pressure_change=False)
    m.fs.Heater_H2 = Heater(property_package=m.fs.props, has_pressure_change=False)
    m.fs.Mixing = Mixer(property_package=m.fs.props, momentum_mixing_type=MomentumMixingType.minimize,
                        inlet_list=["from_o2", "from_h2"])
    m.fs.Mixwater = Mixer(property_package=m.fs.props, momentum_mixing_type=MomentumMixingType.minimize,
                          inlet_list=["from_feed", "from_mix"])

    m.fs.O2_product = Product(property_package=m.fs.props)
    m.fs.H2_product = Product(property_package=m.fs.props)

    # 3. 连接
    m.fs.a01 = Arc(source=m.fs.ely.outlet, destination=m.fs.translator.inlet)
    m.fs.a02 = Arc(source=m.fs.translator.outlet, destination=m.fs.Sep_0.inlet)
    m.fs.a03 = Arc(source=m.fs.Sep_0.o2, destination=m.fs.Sep_O2.inlet)
    m.fs.a04 = Arc(source=m.fs.Sep_0.h2, destination=m.fs.Sep_H2.inlet)
    m.fs.a09 = Arc(source=m.fs.Sep_O2.o2, destination=m.fs.O2_product.inlet)
    m.fs.a10 = Arc(source=m.fs.Sep_H2.h2, destination=m.fs.H2_product.inlet)
    m.fs.a05 = Arc(source=m.fs.Sep_O2.recycle, destination=m.fs.Heater_O2.inlet)
    m.fs.a06 = Arc(source=m.fs.Sep_H2.recycle, destination=m.fs.Heater_H2.inlet)
    m.fs.a07 = Arc(source=m.fs.Heater_O2.outlet, destination=m.fs.Mixing.from_o2)
    m.fs.a08 = Arc(source=m.fs.Heater_H2.outlet, destination=m.fs.Mixing.from_h2)
    m.fs.a11 = Arc(source=m.fs.Mixing.outlet, destination=m.fs.Mixwater.from_mix)
    m.fs.a12 = Arc(source=m.fs.Mixwater.outlet, destination=m.fs.translator1.inlet)
    m.fs.a13 = Arc(source=m.fs.translator1.outlet, destination=m.fs.ely.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m.fs)
    add_scaling(m)
    t = m.fs.time.first()

    # 4. 应用最严格的边界 (消除除以零错误)
    # apply_robust_bounds(m)

    # 5. 参数与固定变量
    N_val, F_val, eta_val, current_val = 326, 96485.33, 0.86, 9135
    prod_H2 = (N_val * current_val * eta_val) / (2 * F_val)
    h2o_consumption = prod_H2

    design_H2O_in = 830
    design_KOH_in = 89
    design_Temp = 343.15
    design_Pres = 16e5
    estimated_split = (0.5 * design_H2O_in - 2 * prod_H2) / (design_H2O_in - h2o_consumption)

    m.fs.ely.outlet.temperature.fix(353.15)
    m.fs.ely.current.fix(current_val)
    m.fs.Heater_O2.outlet.temperature.fix(342.45)
    m.fs.Heater_H2.outlet.temperature.fix(343.75)

    # Translator 约束
    for tr in [m.fs.translator, m.fs.translator1]:
        tr.eq_temp = Constraint(m.fs.time, rule=rule_temp_translation)
        tr.eq_pres = Constraint(m.fs.time, rule=rule_pres_translation)
        tr.eq_flow = ConstraintList()
        for t_idx in m.fs.time:
            for eqn in rule_flow_translation(tr, t_idx):
                tr.eq_flow.add(eqn)

    # 分离设定 (Phase 1)
    m.fs.Sep_0.split_fraction[0, "h2", "Vap", "H2"].fix(1 - 1e-8)
    m.fs.Sep_0.split_fraction[0, "h2", "Vap", "O2"].fix(1e-8)
    m.fs.Sep_0.split_fraction[0, "h2", "Vap", "H2O"].fix(0.67)

    for sep in [m.fs.Sep_O2, m.fs.Sep_H2]:
        sep.split_fraction[0, "recycle", "Liq", "H2O"].fix(1 - 1e-8)
        sep.split_fraction[0, "recycle", "Liq", "KOH"].fix(1 - 1e-8)
        for c in ["H2", "O2", "H2O"]:
            sep.split_fraction[0, "recycle", "Vap", c].fix(1e-8)

    # 耦合约束 (Deactivate)
    def water_split_rule(b, t):
        inlet_water = b.ely.inlet.flow_mol_phase_comp[t, "Liq", "H2O"]
        outlet_water = b.ely.outlet.flow_mol_phase_comp[t, "Liq", "H2O"]
        h2_production = b.ely.outlet.flow_mol_phase_comp[t, "Vap", "H2"]
        target = 0.5 * inlet_water - 2 * h2_production
        return b.Sep_0.split_fraction[t, "h2", "Liq", "H2O"] * outlet_water == target

    m.fs.water_split_constraint = Constraint(m.fs.time, rule=water_split_rule)
    m.fs.water_split_constraint.deactivate()

    def koh_split_rule(b, t):
        return b.Sep_0.split_fraction[t, "h2", "Liq", "KOH"] == b.Sep_0.split_fraction[t, "h2", "Liq", "H2O"]

    m.fs.koh_split_constraint = Constraint(m.fs.time, rule=koh_split_rule)
    m.fs.koh_split_constraint.deactivate()















    # ==================================================================
    # Phase 1: 虚拟流股手动迭代
    # ==================================================================
    print("\n=== Phase 1: Virtual Stream Iteration (Tear a11) ===")

    if hasattr(m.fs.a11, "expanded_block"):
        m.fs.a11.expanded_block.deactivate()
    else:
        raise RuntimeError("Arc expansion failed")

    # 进料固定 (只补水，Phase 1 暂时使用估计值)
    m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, "Liq", "H2O"].fix(h2o_consumption)
    m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, "Liq", "KOH"].fix(1e-12)
    m.fs.Mixwater.from_feed.temperature.fix(design_Temp)
    m.fs.Mixwater.from_feed.pressure.fix(design_Pres)
    for c in ["H2", "O2", "H2O"]: m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, "Vap", c].fix(1e-12)

    # 虚拟流股初始化
    tear_inlet = m.fs.Mixwater.from_mix
    tear_outlet = m.fs.Mixing.outlet

    guess_h2o = design_H2O_in - h2o_consumption
    guess_koh = design_KOH_in

    tear_inlet.unfix()
    tear_inlet.flow_mol_phase_comp[t, "Liq", "H2O"].fix(guess_h2o)
    tear_inlet.flow_mol_phase_comp[t, "Liq", "KOH"].fix(guess_koh)
    tear_inlet.temperature.fix(343.15)
    tear_inlet.pressure.fix(design_Pres)
    for c in ["H2", "O2", "H2O"]: tear_inlet.flow_mol_phase_comp[t, "Vap", c].fix(1e-12)

    m.fs.Sep_0.split_fraction[0, "h2", "Liq", "H2O"].fix(estimated_split)
    m.fs.Sep_0.split_fraction[0, "h2", "Liq", "KOH"].fix(estimated_split)

    solver = get_solver()
    init_level = idaeslog.WARNING

    def initialize_unit(unit):
        try:
            unit.initialize(outlvl=init_level)
        except InitializationError:
            solver.solve(unit, tee=False)
    print("Initializing units...")

    initialize_unit(m.fs.Mixwater)
    propagate_state(m.fs.a12)
    initialize_unit(m.fs.translator1)
    propagate_state(m.fs.a13)
    initialize_unit(m.fs.ely)
    propagate_state(m.fs.a01)
    initialize_unit(m.fs.translator)
    propagate_state(m.fs.a02)
    # initialize_unit(m.fs.Sep_0)
    # propagate_state(m.fs.a03)
    # propagate_state(m.fs.a04)
    # initialize_unit(m.fs.Sep_O2)
    # initialize_unit(m.fs.Sep_H2)
    # propagate_state(m.fs.a05)
    # propagate_state(m.fs.a06)
    # initialize_unit(m.fs.Heater_O2)
    # initialize_unit(m.fs.Heater_H2)
    # propagate_state(m.fs.a07)
    # propagate_state(m.fs.a08)
    # initialize_unit(m.fs.Mixing)

    # 迭代
    print(f"{'Iter':<5} | {'In H2O':<10} | {'Out H2O':<10} | {'Error':<10}")
    damping = 0.5
    for i in range(15):
        try:
            solver.solve(m, tee=False)
        except:
            print("Iteration solve failed.");
            break

        calc_h2o = value(tear_outlet.flow_mol_phase_comp[t, "Liq", "H2O"])
        calc_koh = value(tear_outlet.flow_mol_phase_comp[t, "Liq", "KOH"])

        err = max(abs(calc_h2o - guess_h2o), abs(calc_koh - guess_koh))
        print(f"{i + 1:<5} | {guess_h2o:<10.2f} | {calc_h2o:<10.2f} | {err:<10.4e}")

        if err < 1e-4: break

        guess_h2o += damping * (calc_h2o - guess_h2o)
        guess_koh += damping * (calc_koh - guess_koh)
        tear_inlet.flow_mol_phase_comp[t, "Liq", "H2O"].fix(guess_h2o)
        tear_inlet.flow_mol_phase_comp[t, "Liq", "KOH"].fix(guess_koh)
        tear_inlet.temperature.fix(value(tear_outlet.temperature[t]))

















    # ==================================================================
    # Phase 2: 分步闭环求解 (锚定法)
    # ==================================================================
    print("\n=== Phase 2: Final Closed Loop Solve ===")

    # 1. 连接物理循环
    m.fs.a11.expanded_block.activate()

    # 2. 释放虚拟入口变量 (关键！否则 DOF = -7)
    tear_inlet.unfix()
    print("Unfixed Mixwater.from_mix.")

    # 3. 锚定策略：固定循环流，释放新鲜进料
    #    (让系统自动计算为了维持该循环量需要补充多少水)
    tear_inlet.flow_mol_phase_comp[t, "Liq", "KOH"].fix(guess_koh)
    tear_inlet.flow_mol_phase_comp[t, "Liq", "H2O"].fix(guess_h2o)
    # tear_inlet.pressure.fix(design_Pres)

    m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, "Liq", "KOH"].unfix()
    m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, "Liq", "KOH"].setlb(0)

    m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, "Liq", "H2O"].unfix()
    m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, "Liq", "H2O"].setlb(0)
    # 给一个合理的猜测初值：消耗 + 估算蒸发
    m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, "Liq", "H2O"].value = h2o_consumption + 0.6

    # 4. 激活分离约束
    m.fs.water_split_constraint.activate()
    m.fs.koh_split_constraint.activate()
    m.fs.Sep_0.split_fraction[0, "h2", "Liq", "H2O"].unfix()
    m.fs.Sep_0.split_fraction[0, "h2", "Liq", "KOH"].unfix()

    print(f"Final DOF check: {degrees_of_freedom(m)}")  # 必须为 0

    print("Running Global Solve...")
    solver.options['max_iter'] = 500
    solver.options['tol'] = 1e-6
    # 建议允许求解器容忍微小的计算错误继续尝试
    solver.options['halt_on_ampl_error'] = 'no'

    results = solver.solve(m, tee=True)

    if pyo.check_optimal_termination(results):
        print("\nSUCCESS: Simulation Converged!")
        makeup_h2o = value(m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, 'Liq', 'H2O'])
        print(f"Required Makeup H2O: {makeup_h2o:.4f} mol/s (Evaporation: {makeup_h2o - h2o_consumption:.4f})")
        print(f"Final Split Fraction: {value(m.fs.Sep_0.split_fraction[0, 'h2', 'Liq', 'H2O']):.4f}")
        print(f"Pressure: {value(m.fs.ely.outlet.pressure[t]):.2f} Pa")
    else:
        print("\nWARNING: Solve failed.")
        from idaes.core.util import DiagnosticsToolbox
        dt = DiagnosticsToolbox(m)
        dt.report_structural_issues()
        dt.display_components_with_inconsistent_units()
        dt.display_underconstrained_set()
        dt.display_overconstrained_set()
        dt.display_potential_evaluation_errors()

    analyze_results(m)
    # import os
    # import time
    # if os.path.exists("bihuan.json"):
    #     os.remove("bihuan.json")
    # m.fs.visualize("bihuan")
    # while True:
    #     time.sleep(1)


def initialize_flowsheet_sequentially(m):
    # epsilon = 1e-8
    # # 输入功率
    # u_guess = 1.8
    # current_val = value(m.fs.ely.current[t])
    # N_val = value(m.fs.ely.number_cells)
    # F_val = value(m.fs.ely.F)
    # eta_val = value(m.fs.ely.faraday_efficiency)
    #
    # p_guess = u_guess * current_val * N_val
    # m.fs.ely.voltage[t].set_value(u_guess)
    # m.fs.ely.power[t].set_value(p_guess)
    # if hasattr(m.fs.ely.control_volume, "work"):
    #     m.fs.ely.control_volume.work[t].set_value(p_guess)
    #
    # prod_H2 = (N_val * current_val * eta_val) / (2 * F_val)
    # m.fs.ely.outlet.flow_mol_phase_comp[t, "Vap", "H2"].set_value(prod_H2)
    # m.fs.ely.outlet.flow_mol_phase_comp[t, "Vap", "O2"].set_value(0.5 * prod_H2)
    # F_in_H2O = value(m.fs.water_in.flow_mol_phase_comp[t, "Liq", "H2O"])
    # F_in_KOH = value(m.fs.water_in.flow_mol_phase_comp[t, "Liq", "KOH"])
    # F_vap_H2O_guess = 5
    # m.fs.ely.outlet.flow_mol_phase_comp[t, "Vap", "H2O"].set_value(F_vap_H2O_guess)
    # F_liq_H2O_guess = F_in_H2O - prod_H2 - F_vap_H2O_guess
    # m.fs.ely.outlet.flow_mol_phase_comp[t, "Liq", "H2O"].set_value(max(epsilon, F_liq_H2O_guess))
    # m.fs.ely.outlet.flow_mol_phase_comp[t, "Liq", "KOH"].set_value(F_in_KOH)
    #
    # U_tn = 1.48
    # heat_guess = -(u_guess - U_tn) * current_val * N_val
    # if hasattr(m.fs.ely.control_volume, "heat"):
    #     m.fs.ely.control_volume.heat[t].set_value(heat_guess)
    from idaes.core.util.initialization import propagate_state
    print("Starting sequential initialization...")
    solver = get_solver()
    log = idaeslog.getLogger("Initialization")
    init_level = idaeslog.INFO # 设置为 INFO 以查看详细过程

    # 定义一个辅助函数来进行鲁棒的单元初始化
    def initialize_unit(unit):
        try:
            unit.initialize(outlvl=init_level)
            log.info(f"Successfully initialized {unit.name}")
        except InitializationError:
            # 如果初始化失败，则尝试直接求解该单元块
            log.warning(f"Initialization failed for {unit.name}. Solving block individually.")
            solver.solve(unit, tee=False)

    # 1. Feed (water_in)
    initialize_unit(m.fs.water_in)
    propagate_state(m.fs.a00)
    # 2. 电解槽 (ELY)
    initialize_unit(m.fs.ely)
    propagate_state(m.fs.a01)
    # # 3. Translator
    # initialize_unit(m.fs.translator)
    # propagate_state(m.fs.a02)
    # 4. Sep_0
    # initialize_unit(m.fs.Sep_0)
    # propagate_state(m.fs.a03)
    # propagate_state(m.fs.a04)
    # # 5. Sep_O2 和 Sep_H2
    # initialize_unit(m.fs.Sep_O2)
    # initialize_unit(m.fs.Sep_H2)
    # propagate_state(m.fs.a05)
    # propagate_state(m.fs.a06)
    # 6. Heaters (关键步骤)
    # initialize_unit(m.fs.Heater_O2)
    # initialize_unit(m.fs.Heater_H2)
    # propagate_state(m.fs.a07)
    # propagate_state(m.fs.a08)
    # initialize_unit(m.fs.Mixing)
    print("Flowsheet initialization complete.")




def rule_flow_translation(b, t):
    # 1. 液相 H2O
    yield b.properties_out[t].flow_mol_phase_comp["Liq", "H2O"] == \
        b.properties_in[t].flow_mol_phase_comp["Liq", "H2O"]
    # 2. 液相 KOH
    yield b.properties_out[t].flow_mol_phase_comp["Liq", "KOH"] == \
        b.properties_in[t].flow_mol_phase_comp["Liq", "KOH"]
    # 3. 气相 H2
    yield b.properties_out[t].flow_mol_phase_comp["Vap", "H2"] == \
        b.properties_in[t].flow_mol_phase_comp["Vap", "H2"]
    # 4. 气相 O2
    yield b.properties_out[t].flow_mol_phase_comp["Vap", "O2"] == \
        b.properties_in[t].flow_mol_phase_comp["Vap", "O2"]
    # 5. 气相 H2O
    yield b.properties_out[t].flow_mol_phase_comp["Vap", "H2O"] == \
        b.properties_in[t].flow_mol_phase_comp["Vap", "H2O"]
def rule_temp_translation(b, t):
    return b.properties_out[t].temperature == b.properties_in[t].temperature
def rule_pres_translation(b, t):
    return b.properties_out[t].pressure == b.properties_in[t].pressure
def flow_sheet():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # Property and reaction packages
    m.fs.props = GenericParameterBlock(**props.configuration_VLE)
    m.fs.props0 = GenericParameterBlock(**props.configuration_ELY)
    m.fs.rxn = GenericReactionParameterBlock(property_package=m.fs.props0, **props_reaction.config_dict)

    m.fs.ely = AlkalineElectrolyzer0D(
        property_package=m.fs.props0,
        reaction_package=m.fs.rxn,
        has_heat_transfer=True
        )
    m.fs.translator = Translator(
        inlet_property_package=m.fs.props0,
        outlet_property_package=m.fs.props,
        )

    m.fs.Sep_0 = Separator(
        property_package=m.fs.props,
        outlet_list=["h2", "o2"],
        split_basis=SplittingType.phaseComponentFlow,
        )
    #
    m.fs.Sep_O2 = Separator(
        property_package=m.fs.props,
        outlet_list=["recycle", "o2"],
        split_basis=SplittingType.phaseComponentFlow,
        )
    m.fs.Sep_H2 = Separator(
        property_package=m.fs.props,
        outlet_list=["h2", "recycle"],
        split_basis=SplittingType.phaseComponentFlow,
        )
    #
    #
    m.fs.Heater_O2 = Heater(
        property_package=m.fs.props,
        has_pressure_change=False
        )
    m.fs.Heater_H2 = Heater(
        property_package=m.fs.props,
        has_pressure_change=False
        )


    m.fs.Mixing = Mixer(
        property_package=m.fs.props,
        momentum_mixing_type=MomentumMixingType.minimize,
        inlet_list=["from_o2", "from_h2"],
        )

    m.fs.water_in = Feed(property_package=m.fs.props)
    m.fs.water_out = Product(property_package=m.fs.props)

    m.fs.O2_product = Product(property_package=m.fs.props)
    m.fs.H2_product = Product(property_package=m.fs.props)


    # m.fs.water_feed = Feed(property_package=m.fs.props)

    m.fs.a00 = Arc(source=m.fs.water_in.outlet, destination=m.fs.ely.inlet)

    # m.fs.a0001 = Arc(source=m.fs.water_in.outlet, destination=m.fs.Mixwater.fresh)
    # m.fs.a0002 = Arc(source=m.fs.Mixwater.outlet, destination=m.fs.ely.inlet)
    m.fs.a01 = Arc(source=m.fs.ely.outlet, destination=m.fs.translator.inlet)

    # m.fs.a002 = Arc(source=m.fs.translator.outlet, destination=m.fs.out.inlet)

    m.fs.a02 = Arc(source=m.fs.translator.outlet, destination=m.fs.Sep_0.inlet)


    # m.fs.a0003 = Arc(source=m.fs.Sep_0.o2, destination=m.fs.out1.inlet)
    # m.fs.a0004 = Arc(source=m.fs.Sep_0.h2, destination=m.fs.out2.inlet)

    #
    #
    # # m.fs.a01 = Arc(source=m.fs.ely.outlet, destination=m.fs.Sep_0.inlet)
    #
    m.fs.a03 = Arc(source=m.fs.Sep_0.o2, destination=m.fs.Sep_O2.inlet)
    m.fs.a04 = Arc(source=m.fs.Sep_0.h2, destination=m.fs.Sep_H2.inlet)
    #
    # m.fs.a0010 = Arc(source=m.fs.Sep_O2.o2, destination=m.fs.out1.inlet)
    # m.fs.a0011 = Arc(source=m.fs.Sep_H2.h2, destination=m.fs.out2.inlet)
    # m.fs.a0004 = Arc(source=m.fs.Sep_O2.recycle, destination=m.fs.out3.inlet)
    # m.fs.a0005 = Arc(source=m.fs.Sep_H2.recycle, destination=m.fs.out4.inlet)



    m.fs.a10 = Arc(source=m.fs.Sep_O2.o2, destination=m.fs.O2_product.inlet)
    m.fs.a11 = Arc(source=m.fs.Sep_H2.h2, destination=m.fs.H2_product.inlet)
    m.fs.a05 = Arc(source=m.fs.Sep_O2.recycle, destination=m.fs.Heater_O2.inlet)
    m.fs.a06 = Arc(source=m.fs.Sep_H2.recycle, destination=m.fs.Heater_H2.inlet)

    # m.fs.a0006 = Arc(source=m.fs.Sep_O2.recycle, destination=m.fs.out1.inlet)
    # m.fs.a0007 = Arc(source=m.fs.Sep_H2.recycle, destination=m.fs.out2.inlet)

    # m.fs.a0006 = Arc(source=m.fs.Sep_O2.recycle, destination=m.fs.Mixing.from_o2)
    # m.fs.a0007 = Arc(source=m.fs.Sep_H2.recycle, destination=m.fs.Mixing.from_h2)
    #
    m.fs.a07 = Arc(source=m.fs.Heater_O2.outlet, destination=m.fs.Mixing.from_o2)
    m.fs.a08 = Arc(source=m.fs.Heater_H2.outlet, destination=m.fs.Mixing.from_h2)


    m.fs.a12 = Arc(source=m.fs.Mixing.outlet, destination=m.fs.water_out.inlet)


    # m.fs.a0015 = Arc(source=m.fs.water_in.outlet, destination=m.fs.Mixwater.from_mix)
    # m.fs.a0016 = Arc(source=m.fs.Mixing.outlet, destination=m.fs.Mixwater.recycle)
    TransformationFactory("network.expand_arcs").apply_to(m.fs)
    return m
def set_inputs(m):
    t = m.fs.time.first()
    print("Initial DOF =", degrees_of_freedom(m))
    m.fs.water_in.flow_mol_phase_comp[t, "Liq", "H2O"].fix(830)
    m.fs.water_in.flow_mol_phase_comp[t, "Liq", "KOH"].fix(89)
    m.fs.water_in.flow_mol_phase_comp[t, "Vap", "H2"].fix(1e-8)
    m.fs.water_in.flow_mol_phase_comp[t, "Vap", "O2"].fix(1e-8)
    m.fs.water_in.flow_mol_phase_comp[t, "Vap", "H2O"].fix(1e-8)
    m.fs.water_in.temperature.fix(343.15)
    m.fs.water_in.pressure.fix(16e5)
    print("Initial DOF =", degrees_of_freedom(m))




    m.fs.ely.current.fix(9135)

    m.fs.ely.outlet.temperature.fix(353.15)
    # m.fs.ely.heat_duty.fix(0)







    print("Initial DOF =", degrees_of_freedom(m))

    m.fs.translator.eq_temp = Constraint(m.fs.time, rule=rule_temp_translation)
    m.fs.translator.eq_pres = Constraint(m.fs.time, rule=rule_pres_translation)
    m.fs.translator.eq_flow = ConstraintList()
    for t in m.fs.time:
        for eqn in rule_flow_translation(m.fs.translator, t):
            m.fs.translator.eq_flow.add(eqn)
    print("Initial DOF =", degrees_of_freedom(m))


    m.fs.Sep_0.split_fraction[0, "h2", "Vap", "O2"].fix(1e-8)
    m.fs.Sep_0.split_fraction[0, "h2", "Vap", "H2"].fix(1-1e-8)
    m.fs.Sep_0.split_fraction[0, "h2", "Vap", "H2O"].fix(0.67)
    print("Initial DOF =", degrees_of_freedom(m))
    def water_split_rule(b, t):
        inlet_water = b.ely.inlet.flow_mol_phase_comp[t, "Liq", "H2O"]
        outlet_water = b.ely.outlet.flow_mol_phase_comp[t, "Liq", "H2O"]
        h2_production = b.ely.outlet.flow_mol_phase_comp[t, "Vap", "H2"]
        target_cathode_water = 0.5 * inlet_water - 2 * h2_production
        return b.Sep_0.split_fraction[t, "h2", "Liq", "H2O"] * outlet_water == target_cathode_water
    m.fs.water_split_constraint = Constraint(m.fs.time, rule=water_split_rule)

    def koh_split_rule(b, t):
        return b.Sep_0.split_fraction[t, "h2", "Liq", "KOH"] == \
            b.Sep_0.split_fraction[t, "h2", "Liq", "H2O"]
    m.fs.koh_split_constraint = Constraint(m.fs.time, rule=koh_split_rule)
    print("Initial DOF =", degrees_of_freedom(m))
    #
    m.fs.Sep_O2.split_fraction[0, "recycle", "Vap", "O2"].fix(1e-8)
    m.fs.Sep_O2.split_fraction[0, "recycle", "Vap", "H2"].fix(1e-8)
    m.fs.Sep_O2.split_fraction[0, "recycle", "Vap", "H2O"].fix(1e-8)
    m.fs.Sep_O2.split_fraction[0, "recycle", "Liq", "H2O"].fix(1-1e-8)
    m.fs.Sep_O2.split_fraction[0, "recycle", "Liq", "KOH"].fix(1-1e-8)
    print("Initial DOF =", degrees_of_freedom(m))
    m.fs.Sep_H2.split_fraction[0, "recycle", "Vap", "O2"].fix(1e-8)
    m.fs.Sep_H2.split_fraction[0, "recycle", "Vap", "H2"].fix(1e-8)
    m.fs.Sep_H2.split_fraction[0, "recycle", "Vap", "H2O"].fix(1e-8)
    m.fs.Sep_H2.split_fraction[0, "recycle", "Liq", "H2O"].fix(1-1e-8)
    m.fs.Sep_H2.split_fraction[0, "recycle", "Liq", "KOH"].fix(1-1e-8)
    print("Initial DOF =", degrees_of_freedom(m))

    m.fs.Heater_O2.outlet.temperature.fix(342.45)
    m.fs.Heater_H2.outlet.temperature.fix(343.35)
    print("Initial DOF =", degrees_of_freedom(m))
    return t
def initialize_model(m, t):
    from idaes.core.util.initialization import propagate_state
    print("Starting sequential initialization...")
    solver = get_solver()
    log = idaeslog.getLogger("Initialization")
    init_level = idaeslog.WARNING
    # 定义一个辅助函数
    def initialize_unit(unit):
        try:
            unit.initialize(outlvl=init_level)
            log.info(f"Successfully initialized {unit.name}")
        except InitializationError:
            # 如果初始化失败，则尝试直接求解该单元块
            log.warning(f"Initialization failed for {unit.name}. Solving block individually.")
            solver.solve(unit, tee=False)
    # 1. Feed (water_in)
    initialize_unit(m.fs.water_in)
    propagate_state(m.fs.a00)
    # 2. 电解槽 (ELY)
    initialize_unit(m.fs.ely)
    propagate_state(m.fs.a01)
    # 3. Translator
    initialize_unit(m.fs.translator)
    propagate_state(m.fs.a02)
    # 4. Sep_0
    initialize_unit(m.fs.Sep_0)
    propagate_state(m.fs.a03)
    propagate_state(m.fs.a04)
    # 5. Sep_O2 和 Sep_H2
    initialize_unit(m.fs.Sep_O2)
    initialize_unit(m.fs.Sep_H2)
    propagate_state(m.fs.a05)
    propagate_state(m.fs.a06)
    # 6. Heaters (关键步骤)
    initialize_unit(m.fs.Heater_O2)
    initialize_unit(m.fs.Heater_H2)
    propagate_state(m.fs.a07)
    propagate_state(m.fs.a08)
    initialize_unit(m.fs.Mixing)
def solve(m):
    print(f"DOF : {degrees_of_freedom(m)}")
    solver = pyo.SolverFactory('ipopt')
    solver.options['tol'] = 1e-6
    solver.options['max_iter'] = 300
    print("Starting solve...")
    solver.solve(m, tee=True)
def results(m, t):
    print("\n" + "=" * 40)
    print("Simulation Results (Steady State)")
    print("=" * 40)
    I = pyo.value(m.fs.ely.current[t])
    U = pyo.value(m.fs.ely.voltage[t])
    P = pyo.value(m.fs.ely.power[t])
    T_out = pyo.value(m.fs.ely.outlet.temperature[t])
    H2_prod = pyo.value(m.fs.ely.outlet.flow_mol_phase_comp[t, "Vap", "H2"])
    H2O_in = pyo.value(m.fs.ely.inlet.flow_mol_phase_comp[t, "Liq", "H2O"])
    H2O_out = pyo.value(m.fs.ely.outlet.flow_mol_phase_comp[t, "Liq", "H2O"])
    H2O_consumed = H2O_in - H2O_out
    print(f"Current (I)       : {I:.2f} A")
    print(f"Voltage (U_cell)  : {U:.4f} V")
    print(f"Stack Power       : {P / 1000:.2f} kW")
    print(f"Outlet Temp       : {T_out:.2f} K")
    print("-" * 40)
    print(f"H2 Production     : {H2_prod:.4f} mol/s")
    print(f"Water Consumed    : {H2O_consumed:.4f} mol/s")
    print("-" * 40)
    N_cell = pyo.value(m.fs.ely.number_cells)
    F = pyo.value(m.fs.ely.F)
    z = pyo.value(m.fs.ely.z)
    eta = pyo.value(m.fs.ely.faraday_efficiency)
    theoretical_H2 = (N_cell * I * eta) / (z * F)
    print(f"Calc. H2 Target   : {theoretical_H2:.4f} mol/s")
    H_in = pyo.value(m.fs.ely.control_volume.properties_in[t].enth_mol_phase["Liq"] *
                     m.fs.ely.control_volume.properties_in[t].flow_mol_phase["Liq"])
    H_out = pyo.value(sum(m.fs.ely.control_volume.properties_out[t].enth_mol_phase[p] *
                          m.fs.ely.control_volume.properties_out[t].flow_mol_phase[p] for p in ["Liq", "Vap"]))
    Work_in = pyo.value(m.fs.ely.control_volume.work[t])
    Heat_in = pyo.value(m.fs.ely.control_volume.heat[t])
    print("-" * 40)
    print(f"Energy In (H_in)   : {H_in:.2f} W")
    print(f"Work Input         : {Work_in:.2f} W")
    print(f"Heat Input         : {Heat_in:.2f} W  <-- 检查这一项！")
    print(f"Energy Out (H_out) : {H_out:.2f} W")
    print(f"Balance Error      : {H_in + Work_in + Heat_in - H_out:.2f} W")
    print("-" * 40)
def analyze_results(m):
    t = m.fs.time.first()
    print("\n" + "=" * 80)
    print("RESULTS ANALYSIS REPORT")
    print("=" * 80)

    # --- A. 质量流量统计 (kg/hr) ---
    MW = {'H2': 2.016e-3, 'O2': 31.999e-3, 'H2O': 18.015e-3, 'KOH': 56.105e-3}

    # 定义需要展示的流股
    # 按顺序：电解槽出, 分离器出H2, 分离器出O2, 循环液, 产品H2
    target_arcs = [m.fs.a04, m.fs.a03]

    stream_data = []
    for arc in target_arcs:
        port = arc.source
        name = arc.name

        # 1. 获取基础质量流量
        def get_m(phase, comp):
            try:
                return value(port.flow_mol_phase_comp[t, phase, comp]) * MW[comp] * 3600
            except:
                return 0.0

        m_h2_pure = get_m("Vap", "H2")
        m_o2_pure = get_m("Vap", "O2")
        m_h2o_vap = get_m("Vap", "H2O")
        m_h2o_liq = get_m("Liq", "H2O")
        m_koh = get_m("Liq", "KOH")

        # 2. 归类逻辑
        # 碱液流量
        flow_lye = (m_h2o_liq + m_koh)/3600

        # 气体分配：按 H2/O2 比例分配水蒸气
        total_dry_gas = m_h2_pure + m_o2_pure + 1e-12
        ratio_h2 = m_h2_pure / total_dry_gas
        ratio_o2 = m_o2_pure / total_dry_gas

        flow_h2_stream = (m_h2_pure + m_h2o_vap * ratio_h2)/3600
        flow_o2_stream = (m_o2_pure + m_h2o_vap * ratio_o2)/3600

        stream_data.append({
            "Stream": name,
            "Temp(C)": f"{value(port.temperature[t]) - 273.15:.1f}",
            "Pres(bar)": f"{value(port.pressure[t]) / 1e5:.1f}",
            "H2 Flow (kg/s)": flow_h2_stream,  # 含水蒸气
            "O2 Flow (kg/s)": flow_o2_stream,  # 含水蒸气
            "Lye Flow (kg/s)": flow_lye  # 液体
        })

    print("\n[Mass Flow Statistics]")
    print("Note: Gas flows include associated water vapor.")
    df = pd.DataFrame(stream_data)
    print(df.to_string(index=False, float_format=lambda x: "{:.4f}".format(x)))

    # --- B. 功率与能量平衡 ---
    print("\n[Power & Energy Balance (86% Faraday Eff.)]")

    # 1. 输入电功率
    P_in = value(m.fs.ely.power[t]) / 1000  # kW

    # 2. 有效产出 (H2 LHV)
    LHV_H2_mol = 241.8  # kJ/mol
    n_h2_out = value(m.fs.H2_product.inlet.flow_mol_phase_comp[t, "Vap", "H2"])
    P_out = n_h2_out * LHV_H2_mol  # kW

    # 3. 热损失 (Heat Loss)
    # 电解槽放热 (Q < 0)，取绝对值作为损失
    Q_ely = value(m.fs.ely.control_volume.heat[t]) / 1000
    Loss_Ely = abs(Q_ely) if Q_ely < 0 else 0

    # 换热器冷却负荷
    Q_htr = (value(m.fs.Heater_O2.heat_duty[t]) + value(m.fs.Heater_H2.heat_duty[t])) / 1000
    Loss_Htr = abs(Q_htr) if Q_htr < 0 else 0

    # 4. 显热/其他 (残差)
    Residual = P_in - P_out - Loss_Ely - Loss_Htr

    print(f"{'Item':<30} | {'Value (kW)':<12} | {'Percentage':<10}")
    print("-" * 60)
    print(f"{'Total Electric Power Input':<30} | {P_in:<12.4f} | {'100.00 %':<10}")
    print("-" * 60)
    print(f"{'1. H2 Chemical Energy (LHV)':<30} | {P_out:<12.4f} | {P_out / P_in * 100:<6.2f} %")
    print(f"{'2. Electrolyzer Heat Loss':<30} | {Loss_Ely:<12.4f} | {Loss_Ely / P_in * 100:<6.2f} %")
    print(f"{'3. Heater/Cooler Loss':<30} | {Loss_Htr:<12.4f} | {Loss_Htr / P_in * 100:<6.2f} %")
    print(f"{'4. Sensible Heat / Other':<30} | {Residual:<12.4f} | {Residual / P_in * 100:<6.2f} %")
    print("=" * 80)
    print("Simulation Results (Steady State)")
    print("=" * 40)
    I = pyo.value(m.fs.ely.current[t])
    U = pyo.value(m.fs.ely.voltage[t])
    P = pyo.value(m.fs.ely.power[t])
    print(f"Current (I)       : {I:.2f} A")
    print(f"Voltage (U_cell)  : {U:.4f} V")
    print(f"Stack Power       : {P / 1000:.2f} kW")


def solve_optimization_for_flow(m):
    solver = pyo.SolverFactory('ipopt')
    solver.options['tol'] = 1e-6
    solver.options['max_iter'] = 300
    m.fs.water_in.flow_mol_phase_comp[t, "Liq", "H2O"].unfix()

    m.fs.water_in.flow_mol_phase_comp[t, "Liq", "H2O"].setlb(1.0)

    # 2. 固定（Fix）电解槽出口温度为目标值
    target_temp = 353.15
    # m.fs.ely.outlet.temperature[t].fix(target_temp)
    m.fs.ely.heat_duty[t].fix(0.0)

    # 检查自由度：应该仍然为 0 (Unfix一个, Fix一个)
    print(f"优化模式下的自由度 (DOF): {degrees_of_freedom(m)}")

    if degrees_of_freedom(m) != 0:
        print("警告：自由度不为0，请检查是否有其他变量未正确设置。")

    # ---------------------------------------------------------
    # 第三步：再次求解
    # ---------------------------------------------------------
    print(f"正在求解以匹配目标温度 {target_temp} K ...")
    results_opt = solver.solve(m, tee=True)

    # ---------------------------------------------------------
    # 第四步：输出结果
    # ---------------------------------------------------------
    required_water_flow = value(m.fs.water_in.flow_mol_phase_comp[t, "Liq", "H2O"])
    actual_temp = value(m.fs.ely.outlet.temperature[t])

    print("-" * 30)
    print("计算完成")
    print(f"目标出口温度: {target_temp} K")
    print(f"实际出口温度: {actual_temp:.4f} K")
    print(f"所需的入口水流量 (H2O): {required_water_flow:.4f} mol/s")
    print("-" * 30)

    return required_water_flow


if __name__ == "__main__":
    m = flow_sheet()
    t = set_inputs(m)
    initialize_model(m, t)

    solve(m)
    # solve_optimization_for_flow(m)


    # results(m, t)
    # analyze_results(m)

    import os
    import time
    if os.path.exists("kaihuan.json"):
        os.remove("kaihuan.json")
    m.fs.visualize("kaihuan")
    while True:
        time.sleep(1)

    # from idaes.core.util import DiagnosticsToolbox
    # dt = DiagnosticsToolbox(m)
    # dt.report_structural_issues()
    # dt.display_components_with_inconsistent_units()
    # dt.display_underconstrained_set()
    # dt.display_overconstrained_set()
    # dt.display_potential_evaluation_errors()






