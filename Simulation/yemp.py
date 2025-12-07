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
import matplotlib
matplotlib.use('TkAgg')  # 或者 'Qt5Agg', 'Agg'
import matplotlib.pyplot as plt
import numpy as np
from idaes.core.solvers import get_solver

import props as props
import props_reaction as props_reaction
from idaes.core.util.initialization import propagate_state
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
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor
from idaes.core import StateBlock
from idaes.core import UnitModelBlockData
import logging
import idaes.logger as idaeslog
from idaes.core import StateBlock, ControlVolume0DBlock
from idaes.core.util import scaling as iscale
import pandas as pd
from pyomo.opt import SolverStatus, TerminationCondition, check_optimal_termination


def create_flowsheet():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = GenericParameterBlock(**props.configuration_VLE)
    m.fs.props0 = GenericParameterBlock(**props.configuration_ELY)
    m.fs.rxn = GenericReactionParameterBlock(property_package=m.fs.props0, **props_reaction.config_dict)

    m.fs.ely = AlkalineElectrolyzer0D(property_package=m.fs.props0, reaction_package=m.fs.rxn,         has_heat_transfer=True)
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
    for v in m.component_data_objects(Var, descend_into=True):
        if "flow_mol_phase_comp" in v.name: set_scaling_factor(v, 1e-2)
    calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    t = m.fs.time.first()

    m.params = {}
    m.params['N_val'] = 326
    m.params['F_val'] = 96485.33
    m.params['eta_val'] = 0.86
    m.params['current_val'] = 9135

    prod_H2 = (m.params['N_val'] * m.params['current_val'] * m.params['eta_val']) / (2 * m.params['F_val'])
    m.params['h2o_consumption'] = prod_H2

    design_H2O_in = 830
    m.params['estimated_split'] = (0.5 * design_H2O_in - 2 * prod_H2) / (design_H2O_in - m.params['h2o_consumption'])

    # m.fs.ely.outlet.temperature.fix(353.15)
    m.fs.ely.heat_duty.fix(0)
    m.fs.ely.current.fix(m.params['current_val'])
    m.fs.Heater_O2.outlet.temperature.fix(342.45)
    m.fs.Heater_H2.outlet.temperature.fix(343.75)

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

    for tr in [m.fs.translator, m.fs.translator1]:
        tr.eq_temp = Constraint(m.fs.time, rule=rule_temp_translation)
        tr.eq_pres = Constraint(m.fs.time, rule=rule_pres_translation)
        tr.eq_flow = ConstraintList()
        for t_idx in m.fs.time:
            for eqn in rule_flow_translation(tr, t_idx):
                tr.eq_flow.add(eqn)

    m.fs.Sep_0.split_fraction[0, "h2", "Vap", "H2"].fix(1 - 1e-8)
    m.fs.Sep_0.split_fraction[0, "h2", "Vap", "O2"].fix(1e-8)
    m.fs.Sep_0.split_fraction[0, "h2", "Vap", "H2O"].fix(0.67)  # 初始猜测

    for sep in [m.fs.Sep_O2, m.fs.Sep_H2]:
        sep.split_fraction[0, "recycle", "Liq", "H2O"].fix(1 - 1e-8)
        sep.split_fraction[0, "recycle", "Liq", "KOH"].fix(1 - 1e-8)
        for c in ["H2", "O2", "H2O"]:
            sep.split_fraction[0, "recycle", "Vap", c].fix(1e-8)

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


def initialize_system(m, design_H2O_in=830, design_KOH_in=89, design_Temp=343.15, design_Pres=16e5, damping=0.5, max_iter=30, force_run=False):
    print("\n=== 第一步: 虚拟流股撕裂 (a11) ===")
    t = m.fs.time.first()
    solver = get_solver()

    # 1. 撕裂流股
    if hasattr(m.fs.a11, "expanded_block"):
        m.fs.a11.expanded_block.deactivate()
    else:
        raise RuntimeError("Arc expansion failed")

    # 2. 进料固定 (只补水，Phase 1 暂时使用估计值)
    h2o_consumption = m.params['h2o_consumption']

    m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, "Liq", "H2O"].fix(h2o_consumption)
    m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, "Liq", "KOH"].fix(1e-12)
    m.fs.Mixwater.from_feed.temperature.fix(design_Temp)
    m.fs.Mixwater.from_feed.pressure.fix(design_Pres)
    for c in ["H2", "O2", "H2O"]:
        m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, "Vap", c].fix(1e-12)

    # 3. 虚拟流股初值设置
    tear_inlet = m.fs.Mixwater.from_mix
    tear_outlet = m.fs.Mixing.outlet

    guess_h2o = design_H2O_in - h2o_consumption
    guess_koh = design_KOH_in

    tear_inlet.unfix()
    tear_inlet.flow_mol_phase_comp[t, "Liq", "H2O"].fix(guess_h2o)
    tear_inlet.flow_mol_phase_comp[t, "Liq", "KOH"].fix(guess_koh)
    tear_inlet.temperature.fix(343.15)
    tear_inlet.pressure.fix(design_Pres)
    for c in ["H2", "O2", "H2O"]:
        tear_inlet.flow_mol_phase_comp[t, "Vap", c].fix(1e-12)

    # 设置分离器初值
    m.fs.Sep_0.split_fraction[0, "h2", "Liq", "H2O"].fix(m.params['estimated_split'])
    m.fs.Sep_0.split_fraction[0, "h2", "Liq", "KOH"].fix(m.params['estimated_split'])


    print("Initializing units...")
    init_level = idaeslog.ERROR
    def initialize_unit(unit):
        try:
            unit.initialize(outlvl=init_level)
        except InitializationError:
            solver.solve(unit, tee=False)
    initialize_unit(m.fs.Mixwater)
    propagate_state(m.fs.a12)
    initialize_unit(m.fs.translator1)
    propagate_state(m.fs.a13)
    initialize_unit(m.fs.ely)
    propagate_state(m.fs.a01)

    initialize_unit(m.fs.translator)
    propagate_state(m.fs.a02)  # 将状态传给 Sep_0

    initialize_unit(m.fs.Sep_0)
    propagate_state(m.fs.a03)  # -> Sep_O2
    propagate_state(m.fs.a04)  # -> Sep_H2

    initialize_unit(m.fs.Sep_O2)
    propagate_state(m.fs.a09)  # -> O2_product
    propagate_state(m.fs.a05)  # -> Heater_O2

    initialize_unit(m.fs.O2_product)

    initialize_unit(m.fs.Heater_O2)
    propagate_state(m.fs.a07)  # -> Mixing (from_o2)

    initialize_unit(m.fs.Sep_H2)
    propagate_state(m.fs.a10)  # -> H2_product
    propagate_state(m.fs.a06)  # -> Heater_H2

    initialize_unit(m.fs.H2_product)

    initialize_unit(m.fs.Heater_H2)
    propagate_state(m.fs.a08)  # -> Mixing (from_h2)



    print(f"{'Iter':<5} | {'In H2O':<10} | {'Out H2O':<10} | {'Error':<10}")
    damping = damping
    for i in range(max_iter, ):
        try:
            solver.solve(m, tee=False)
        except:
            print("Iteration solve failed.")
            break
        calc_h2o = value(tear_outlet.flow_mol_phase_comp[t, "Liq", "H2O"])
        calc_koh = value(tear_outlet.flow_mol_phase_comp[t, "Liq", "KOH"])
        err = max(abs(calc_h2o - guess_h2o), abs(calc_koh - guess_koh))
        print(f"{i + 1:<5} | {guess_h2o:<10.2f} | {calc_h2o:<10.2f} | {err:<10.4e}")
        if not force_run and err < 1e-4:
            break

        guess_h2o += damping * (calc_h2o - guess_h2o)
        guess_koh += damping * (calc_koh - guess_koh)
        tear_inlet.flow_mol_phase_comp[t, "Liq", "H2O"].fix(guess_h2o)
        tear_inlet.flow_mol_phase_comp[t, "Liq", "KOH"].fix(guess_koh)
        tear_inlet.temperature.fix(value(tear_outlet.temperature[t]))

    m.params['converged_h2o'] = guess_h2o
    m.params['converged_koh'] = guess_koh


def solve_closed_loop(m):
    print("\n=== 第二步: 循环求解 ===")
    t = m.fs.time.first()
    solver = get_solver()

    tear_inlet = m.fs.Mixwater.from_mix

    m.fs.a11.expanded_block.activate()

    tear_inlet.unfix()
    print("Unfixed Mixwater.from_mix.")

    # 3. 锚定策略：固定循环流，释放新鲜进料
    # 使用 Phase 1 收敛得到的值
    tear_inlet.flow_mol_phase_comp[t, "Liq", "KOH"].fix(m.params['converged_koh'])
    tear_inlet.flow_mol_phase_comp[t, "Liq", "H2O"].fix(m.params['converged_h2o'])

    m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, "Liq", "KOH"].unfix()
    m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, "Liq", "KOH"].setlb(0)

    m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, "Liq", "H2O"].unfix()
    m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, "Liq", "H2O"].setlb(0)

    # 猜测初值
    m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, "Liq", "H2O"].value = m.params['h2o_consumption'] + 0.6

    # 4. 激活分离约束
    m.fs.water_split_constraint.activate()
    m.fs.koh_split_constraint.activate()
    m.fs.Sep_0.split_fraction[0, "h2", "Liq", "H2O"].unfix()
    m.fs.Sep_0.split_fraction[0, "h2", "Liq", "KOH"].unfix()

    print(f"Final DOF check: {degrees_of_freedom(m)}")

    print("Running Global Solve...")
    solver.options['max_iter'] = 500
    solver.options['tol'] = 1e-6
    solver.options['halt_on_ampl_error'] = 'no'

    results = solver.solve(m, tee=False)

    if pyo.check_optimal_termination(results):
        print("\nSUCCESS: Simulation Converged!")
        makeup_h2o = value(m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, 'Liq', 'H2O'])
        print(f"Required Makeup H2O: {makeup_h2o:.4f} mol/s")
    else:
        print("\nWARNING: Solve failed.")
    return results


def analyze_results(m):
    t = m.fs.time.first()
    print("\n" + "=" * 80)
    print("RESULTS ANALYSIS REPORT")
    print("=" * 80)

    # --- A. 质量流量统计 (kg/hr) ---
    MW = {'H2': 2.016e-3, 'O2': 31.999e-3, 'H2O': 18.015e-3, 'KOH': 56.105e-3}

    # 定义需要展示的流股
    # 按顺序：电解槽出, 分离器出H2, 分离器出O2, 循环液, 产品H2
    target_arcs = [m.fs.a01, m.fs.a04, m.fs.a03, m.fs.a11, m.fs.a10, m.fs.a09]

    stream_data = []
    for arc in target_arcs:
        port = arc.source
        name = arc.name

    #     # 1. 获取基础质量流量
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
    LHV_H2_mol = 285.8  # kJ/mol
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

    print(f"H2 Production     : {n_h2_out:.4f} mol/s")
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


def run():
    m = create_flowsheet()
    set_operating_conditions(m)
    initialize_system(m)
    solve_closed_loop(m)
    analyze_results(m)
    import os
    import time
    if os.path.exists("xunhuan.json"):
        os.remove("xunhuan.json")
    m.fs.visualize("xunhuan")
    while True:
        time.sleep(1)

def diedai_canshu():
    damping_list = [0.3, 0.4, 0.5, 0.6, 0.7]
    iter_list = list(range(10, 46, 5))
    # damping_list = [0.3]
    # iter_list = [10]
    results_data = {d: [] for d in damping_list}

    print("==================================================")
    print("Starting Parameter Sensitivity Analysis")
    print("==================================================")
    print(f"{'Damping':<10} | {'Iter':<10} | {'Makeup H2O (kg/s)':<20} | {'Status'}")
    print("-" * 60)

    # 双重循环遍历参数
    for damp in damping_list:
        for n_iter in iter_list:

            m = create_flowsheet()
            set_operating_conditions(m)
            result_kg_s = np.nan
            status = "Failed"

            try:
                initialize_system(m, damping=damp, max_iter=n_iter, force_run=True)
                solve_results = solve_closed_loop(m)
                if check_optimal_termination(solve_results):
                    t = m.fs.time.first()
                    val_mol_s = value(m.fs.Mixwater.from_feed.flow_mol_phase_comp[t, 'Liq', 'H2O'])
                    result_kg_s = val_mol_s * 0.018015
                    status = "Optimal"
                else:
                    status = "No optimal solution"
            except Exception as e:
                status = f"Error: {str(e)[:20]}..."

            print(f"{damp:<10.1f} | {n_iter:<10d} | {result_kg_s:<20.6f} | {status}")
            results_data[damp].append(result_kg_s)

    print("\nAnalysis Complete. Generating Plot...")

    plt.figure(figsize=(12, 7))
    markers = ['o', 's', '^', 'D', 'v']

    for i, damp in enumerate(damping_list):
        plt.plot(iter_list, results_data[damp],
                 marker=markers[i % len(markers)],
                 markersize=8,
                 linewidth=2,
                 label=f'Damping = {damp}')

    plt.xlabel('Initialization Iterations', fontsize=12, fontweight='bold')
    plt.ylabel('Final Makeup Water Flow (kg/s)', fontsize=12, fontweight='bold')
    plt.title('Impact of Initialization Parameters on Final Solution', fontsize=14)
    plt.grid(True, which='major', linestyle='-', alpha=0.6)
    plt.grid(True, which='minor', linestyle='--', alpha=0.3)
    plt.minorticks_on()
    plt.legend(fontsize=11)

    plt.tight_layout()

    plt.savefig('sensitivity_analysis_result.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    idaeslog.getLogger("idaes").setLevel(idaeslog.ERROR)
    logging.getLogger('pyomo').setLevel(logging.ERROR)
    logging.getLogger("idaes.core.util.scaling").setLevel(logging.ERROR)
    run()

    # import os
    # import time
    # if os.path.exists("xunhuan.json"):
    #     os.remove("xunhuan.json")
    # m.fs.visualize("xunhuan")
    # while True:
    #     time.sleep(1)





