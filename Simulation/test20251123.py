import pyomo.environ as pyo
from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock
)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock
)

# 导入你之前建立的模块
import props
import props_reaction
from Electronchemical0D import AlkalineElectrolyzer0D
from pyomo.environ import value
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor
from idaes.core.util.model_statistics import degrees_of_freedom
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
def main000():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = GenericParameterBlock(**props.configuration_ELY)
    m.fs.rxn = GenericReactionParameterBlock(property_package=m.fs.props, **props_reaction.config_dict)
    m.fs.ely = AlkalineElectrolyzer0D(property_package=m.fs.props, reaction_package=m.fs.rxn,         has_heat_transfer=True)
    print(f"Initial DOF: {degrees_of_freedom(m)}")
    t = m.fs.time.first()
    # --- 1. 设定操作参数 ---
    m.fs.ely.inlet.temperature.fix(343.15)
    m.fs.ely.inlet.pressure.fix(16e5)
    print(f"Initial DOF: {degrees_of_freedom(m)}")
    m.fs.ely.inlet.flow_mol_phase_comp[t, "Liq", "H2O"].fix(830)
    m.fs.ely.inlet.flow_mol_phase_comp[t, "Liq", "KOH"].fix(89)
    m.fs.ely.inlet.flow_mol_phase_comp[t, "Vap", "H2"].fix(1e-6)
    m.fs.ely.inlet.flow_mol_phase_comp[t, "Vap", "O2"].fix(1e-6)
    m.fs.ely.inlet.flow_mol_phase_comp[t, "Vap", "H2O"].fix(1e-6)
    m.fs.ely.current.fix(9135)
    m.fs.ely.outlet.temperature.fix(353.15)
    print(f"Initial DOF: {degrees_of_freedom(m)}")

    import idaes.logger as idaeslog
    init_level = idaeslog.ERROR
    def initialize_unit(unit):
        unit.initialize(outlvl=init_level)
    initialize_unit(m.fs.ely)

    solver = pyo.SolverFactory('ipopt')
    solver.options['tol'] = 1e-5
    solver.options['max_iter'] = 100
    solver.options['nlp_scaling_method'] = 'gradient-based'
    print("Starting solve...")
    results = solver.solve(m, tee=True)
    U = pyo.value(m.fs.ely.voltage[t])
    H2_prod = pyo.value(m.fs.ely.outlet.flow_mol_phase_comp[t, "Vap", "H2"])
    print(f"Voltage (U_cell)  : {U:.4f} V")
    print(f"H2 Production     : {H2_prod:.4f} mol/s")





    # #
    # m.fs.ely.control_volume.initialize()
    # # 输入功率
    # u_guess = 1.8
    # current_val = value(m.fs.ely.current[t])
    # p_guess = u_guess * current_val * value(m.fs.ely.number_cells)  # W
    # m.fs.ely.voltage[t].set_value(u_guess)
    # m.fs.ely.power[t].set_value(p_guess)
    # if hasattr(m.fs.ely.control_volume, "work"):
    #     m.fs.ely.control_volume.work[t].set_value(p_guess)
    # # H2产量
    # F = value(m.fs.ely.F)
    # N = value(m.fs.ely.number_cells)
    # eta = value(m.fs.ely.faraday_efficiency)
    # prod_H2 = (N * current_val * eta) / (2 * F)
    # m.fs.ely.outlet.flow_mol_phase_comp[t, "Vap", "H2"].set_value(prod_H2)
    # m.fs.ely.outlet.flow_mol_phase_comp[t, "Vap", "O2"].set_value(0.5 * prod_H2)
    # m.fs.ely.outlet.flow_mol_phase_comp[t, "Liq", "H2O"].set_value(100 - prod_H2)
    # m.fs.ely.outlet.flow_mol_phase_comp[t, "Liq", "KOH"].set_value(1)
    #
    # # 放热
    # heat_guess = (1.48 - u_guess) * current_val * N
    # if hasattr(m.fs.ely.control_volume, "heat"):
    #     m.fs.ely.control_volume.heat[t].set_value(heat_guess)


    # 6. 结果输出与验证
    # print("\n" + "=" * 40)
    # print("Simulation Results (Steady State)")
    # print("=" * 40)
    #
    # # 获取关键变量值
    # I = pyo.value(m.fs.ely.current[t])
    # U = pyo.value(m.fs.ely.voltage[t])
    # P = pyo.value(m.fs.ely.power[t])
    # T_out = pyo.value(m.fs.ely.outlet.temperature[t])
    #
    # H2_prod = pyo.value(m.fs.ely.outlet.flow_mol_phase_comp[t, "Vap", "H2"])
    # H2O_in = pyo.value(m.fs.ely.inlet.flow_mol_phase_comp[t, "Liq", "H2O"])
    # H2O_out = pyo.value(m.fs.ely.outlet.flow_mol_phase_comp[t, "Liq", "H2O"])
    # H2O_consumed = H2O_in - H2O_out
    #
    # print(f"Current (I)       : {I:.2f} A")
    # print(f"Voltage (U_cell)  : {U:.4f} V")
    # print(f"Stack Power       : {P / 1000:.2f} kW")
    # print(f"Outlet Temp       : {T_out:.2f} K")
    # print("-" * 40)
    # print(f"H2 Production     : {H2_prod:.4f} mol/s")
    # print(f"Water Consumed    : {H2O_consumed:.4f} mol/s")
    # print("-" * 40)
    #
    # N_cell = pyo.value(m.fs.ely.number_cells)
    # F = pyo.value(m.fs.ely.F)
    # z = pyo.value(m.fs.ely.z)
    # eta = pyo.value(m.fs.ely.faraday_efficiency)
    #
    # theoretical_H2 = (N_cell * I * eta) / (z * F)
    # print(f"Calc. H2 Target   : {theoretical_H2:.4f} mol/s")
    # # 打印能量项
    # H_in = pyo.value(m.fs.ely.control_volume.properties_in[t].enth_mol_phase["Liq"] *
    #                  m.fs.ely.control_volume.properties_in[t].flow_mol_phase["Liq"])
    # H_out = pyo.value(sum(m.fs.ely.control_volume.properties_out[t].enth_mol_phase[p] *
    #                       m.fs.ely.control_volume.properties_out[t].flow_mol_phase[p] for p in ["Liq", "Vap"]))
    # Work_in = pyo.value(m.fs.ely.control_volume.work[t])
    # Heat_in = pyo.value(m.fs.ely.control_volume.heat[t])
    #
    # print("-" * 40)
    # print(f"Energy In (H_in)   : {H_in:.2f} W")
    # print(f"Work Input         : {Work_in:.2f} W")
    # print(f"Heat Input         : {Heat_in:.2f} W  <-- 检查这一项！")
    # print(f"Energy Out (H_out) : {H_out:.2f} W")
    # print(f"Balance Error      : {H_in + Work_in + Heat_in - H_out:.2f} W")
    # print("-" * 40)
    # if abs(H2_prod - theoretical_H2) < 1e-3:
    #     print(" Test PASSED: Mass balance matches Faraday's Law.")
    #     import os
    #     import time
    #     if os.path.exists("dianjiecao.json"):
    #         os.remove("dianjiecao.json")
    #     m.fs.visualize("dianjiecao")
    #     while True:
    #         time.sleep(1)
    # else:
    #     print(" Test FAILED: Production rate mismatch.")

    # from idaes.core.util import DiagnosticsToolbox
    # dt = DiagnosticsToolbox(m)
    # dt.report_structural_issues()
    # dt.display_components_with_inconsistent_units()
    # dt.display_underconstrained_set()
    # dt.display_overconstrained_set()
    # dt.display_potential_evaluation_errors()
def main():
    # 存储结果的列表
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']  # 设置中文字体
    plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号
    current_densities = []
    voltages = []

    # 电流密度范围：0到0.35 A/cm²
    # 转换为A/m²：1 A/cm² = 10000 A/m²
    current_density_range = np.linspace(0.001, 0.35, 35)  # 从0.001开始，避免0电流的数值问题

    for j in current_density_range:
        print(f"\n=== 计算电流密度: {j:.3f} A/cm2 ===")

        # 创建新模型实例
        m = pyo.ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.props = GenericParameterBlock(**props.configuration_ELY)
        m.fs.rxn = GenericReactionParameterBlock(property_package=m.fs.props,
                                                 **props_reaction.config_dict)
        m.fs.ely = AlkalineElectrolyzer0D(
            property_package=m.fs.props,
            reaction_package=m.fs.rxn,
        has_heat_transfer=True
        )

        t = m.fs.time.first()

        # --- 1. 设定操作参数 ---
        m.fs.ely.inlet.temperature.fix(343.15)  # 70°C
        m.fs.ely.inlet.pressure.fix(16e5)  # 16 bar

        # 设定进料
        m.fs.ely.inlet.flow_mol_phase_comp[t, "Liq", "H2O"].fix(100)
        m.fs.ely.inlet.flow_mol_phase_comp[t, "Liq", "KOH"].fix(1)
        m.fs.ely.inlet.flow_mol_phase_comp[t, "Vap", "H2"].fix(1e-6)
        m.fs.ely.inlet.flow_mol_phase_comp[t, "Vap", "O2"].fix(1e-6)
        m.fs.ely.inlet.flow_mol_phase_comp[t, "Vap", "H2O"].fix(1e-6)

        # --- 2. 根据电流密度计算并设定电流 ---
        # 电流密度 j (A/cm²) -> 电流 I (A)
        # 面积: 2.66 m² = 26600 cm²
        electrode_area_cm2 = 2.66 * 10000  # m² 转 cm²
        current = j * electrode_area_cm2  # I = j × A
        print(f"电流密度: {j:.3f} A/cm2")
        print(f"对应电流: {current:.1f} A")

        m.fs.ely.current.fix(current)
        m.fs.ely.outlet.temperature.fix(353.15)  # 80°C

        # --- 3. 初始化并求解 ---
        import idaes.logger as idaeslog
        init_level = idaeslog.ERROR

        def initialize_unit(unit):
            unit.initialize(outlvl=init_level)

        initialize_unit(m.fs.ely)

        solver = pyo.SolverFactory('ipopt')
        solver.options['tol'] = 1e-5
        solver.options['max_iter'] = 100
        solver.options['nlp_scaling_method'] = 'gradient-based'

        print(f"求解电流密度 {j:.3f} A/cm²...")
        results = solver.solve(m, tee=False)  # 设为False避免输出过多

        if results.solver.termination_condition == pyo.TerminationCondition.optimal:
            voltage = pyo.value(m.fs.ely.voltage[t])
            print(f"电压: {voltage:.3f} V")

            # 保存结果
            current_densities.append(j)
            voltages.append(voltage)
        else:
            print(f"警告: 电流密度 {j:.3f} A/cm² 求解失败")
            # 可以记录NaN或跳过

    # --- 4. 绘制V-j曲线 ---
    plt.figure(figsize=(10, 6))
    plt.plot(current_densities, voltages, 'bo-', linewidth=2, markersize=8,
             label='V-j 特性曲线')

    # 添加理论分解电压线（约1.23V）
    plt.axhline(y=1.23, color='r', linestyle='--', alpha=0.5,
                label='理论分解电压 (1.23V)')

    plt.xlabel('电流密度 (A/cm2)', fontsize=12)
    plt.ylabel('槽电压 (V)', fontsize=12)
    plt.title('碱性电解槽 V-j 特性曲线\n温度: 70-80°C, 压力: 16 bar', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)

    # 设置坐标轴范围
    plt.xlim(0, max(current_densities) * 1.05)
    plt.ylim(min(voltages) * 0.95, max(voltages) * 1.05)

    # 添加数据点标记
    for i, (j, v) in enumerate(zip(current_densities, voltages)):
        if i % 5 == 0:  # 每5个点标记一次
            plt.annotate(f'{v:.2f}V', xy=(j, v), xytext=(5, 5),
                         textcoords='offset points', fontsize=9)

    plt.tight_layout()

    # 保存图像
    plt.savefig('voltage_current_density_curve.png', dpi=300, bbox_inches='tight')
    plt.show()

    # --- 5. 打印结果表格 ---
    print("\n" + "=" * 50)
    print("计算结果汇总:")
    print("=" * 50)
    print(f"{'电流密度 (A/cm2)':<20} {'电流 (A)':<15} {'电压 (V)':<10}")
    print("-" * 50)

    for j, v in zip(current_densities, voltages):
        current = j * electrode_area_cm2
        print(f"{j:<20.3f} {current:<15.1f} {v:<10.3f}")

    return current_densities, voltages

def load_literature_data(csv_path):
    """
    从CSV加载文献数据（已自动按横坐标排序）
    csv_path: CSV文件路径，第一列横坐标，第二列纵坐标
    """
    df = pd.read_csv(csv_path, header=None)  # 无表头
    literature_data = {
        'current_density': df[0].tolist(),  # 横坐标：电流密度
        'voltage': df[1].tolist()  # 纵坐标：电压
    }
    # 重要：按横坐标排序（解决顺序问题）
    sort_idx = np.argsort(literature_data['current_density'])
    literature_data['current_density'] = [literature_data['current_density'][i] for i in sort_idx]
    literature_data['voltage'] = [literature_data['voltage'][i] for i in sort_idx]

    print(f"从 {csv_path} 加载了 {len(literature_data['current_density'])} 个数据点（已按横坐标排序）")
    return literature_data
def plot_comparison_only(literature_csv_path, your_current, your_voltage,
                         electrode_area=2.66, save_path='comparison.png'):
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']  # 设置中文字体
    plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

    literature = load_literature_data(literature_csv_path)
    lit_x = np.array(literature['current_density'])
    lit_y = np.array(literature['voltage'])


    your_x = np.array(your_current)
    your_y = np.array(your_voltage)

    # 重要：对您的数据也按横坐标排序（确保顺序一致）
    sort_idx_own = np.argsort(your_x)
    your_x = your_x[sort_idx_own]
    your_y = your_y[sort_idx_own]

    # 3. 绘图
    plt.figure(figsize=(10, 6))

    # 文献数据 - 红色实线
    plt.plot(lit_x, lit_y,
             color='red',
             marker='o',
             linestyle='-',
             linewidth=2,
             markersize=8,
             markerfacecolor='white',
             markeredgewidth=2,
             label='文献数据',
             alpha=0.8)

    # 您的数据 - 蓝色虚线
    plt.plot(your_x, your_y,
             color='blue',
             marker='s',
             linestyle='--',
             linewidth=2,
             markersize=8,
             markerfacecolor='white',
             markeredgewidth=2,
             label='计算结果',
             alpha=0.8)

    # 4. 图表设置
    plt.xlabel('电流密度 (A/cm2)', fontsize=12)
    plt.ylabel('槽电压 (V)', fontsize=12)
    plt.title('V-j 特性曲线对比', fontsize=14)
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.legend(loc='best', fontsize=11)

    # 5. 自动调整坐标轴
    all_x = np.concatenate([lit_x, your_x])
    plt.xlim(0, max(all_x) * 1.05)

    all_y = np.concatenate([lit_y, your_y])
    y_min, y_max = min(all_y), max(all_y)
    margin = (y_max - y_min) * 0.1
    plt.ylim(y_min - margin, y_max + margin)

    # 6. 保存和显示
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

    print(f"对比图已保存: {save_path}")
    print(f"文献数据点: {len(lit_x)} 个（已排序）")
    print(f"您的数据点: {len(your_x)} 个（已排序）")
def V_currentdensity(current_densities, voltages):
    your_data = {
        'current': list(current_densities),  # 转为列表
        'voltage': list(voltages)
    }
    literature_csv = "literature.csv"
    plot_comparison_only(
        literature_csv_path=literature_csv,
        your_current=your_data['current'],
        your_voltage=your_data['voltage'],
        save_path="result.png"
    )





def calculate_h2_production(current_range):
    h2_production = []
    for current in current_range:
        faraday = 96485.33  # C/mol
        h2_mol_per_sec = (163*current*0.86 / (2 * faraday))  # mol/s
        h2_nmc3_per_hour = h2_mol_per_sec * 3600 * 0.022414
        h2_production.append(h2_nmc3_per_hour)
    return h2_production
def load_h2_literature_data(csv_path):
    df = pd.read_csv(csv_path, header=None)  # 无表头
    literature_data = {
        'current': df[0].tolist(),  # 横坐标：电流(A)
        'h2_production': df[1].tolist()  # 纵坐标：氢气产量(Nm³/h)
    }
    sort_idx = np.argsort(literature_data['current'])
    literature_data['current'] = [literature_data['current'][i] for i in sort_idx]
    literature_data['h2_production'] = [literature_data['h2_production'][i] for i in sort_idx]
    print(f"从 {csv_path} 加载了 {len(literature_data['current'])} 个数据点（已按电流排序）")
    return literature_data
def plot_h2_production_comparison(literature_csv_path,calculated_current,calculated_h2,save_path='h2_production_comparison.png'):
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']  # 设置中文字体
    plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号
    literature = load_h2_literature_data(literature_csv_path)
    lit_x = np.array(literature['current'])
    lit_y = np.array(literature['h2_production'])
    calc_x = np.array(calculated_current)
    calc_y = np.array(calculated_h2)
    plt.figure(figsize=(10, 6))
    # 文献数据 - 红色实线
    plt.plot(lit_x, lit_y,
             color='red',
             marker='o',
             linestyle='-',
             linewidth=2,
             markersize=8,
             markerfacecolor='white',
             markeredgewidth=2,
             label='文献数据',
             alpha=0.8)
    # 计算数据 - 蓝色虚线
    plt.plot(calc_x, calc_y,
             color='blue',
             marker='s',
             linestyle='--',
             linewidth=2,
             markersize=8,
             markerfacecolor='white',
             markeredgewidth=2,
             label='计算结果',
             alpha=0.8)
    # 设置图表
    plt.xlabel('电流 (A)', fontsize=12)
    plt.ylabel('氢气产量 (Nm3/h)', fontsize=12)
    plt.title('氢气产量-电流特性曲线对比', fontsize=14)
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.legend(loc='best', fontsize=11)

    # 设置坐标轴范围
    plt.xlim(5500, 9500)  # 电流范围
    plt.ylim(0, 700)  # 氢气产量范围

    # 保存和显示
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

    print(f"氢气产量对比图已保存: {save_path}")
    print(f"文献数据点: {len(lit_x)} 个（已排序）")
    print(f"计算数据点: {len(calc_x)} 个")
def H2_Current():
    current_range = np.arange(5500, 9501, 500)
    calculated_h2 = calculate_h2_production(current_range)
    literature_csv = "literature_h2.csv"
    plot_h2_production_comparison(
        literature_csv_path=literature_csv,
        calculated_current=current_range,
        calculated_h2=calculated_h2,
        save_path="h2_production_comparison.png"
    )
    print("\n计算结果 (电流 vs 氢气产量):")
    for i, current in enumerate(current_range):
        print(f"电流: {current} A → 氢气产量: {calculated_h2[i]:.2f} Nm3/h")
if __name__ == "__main__":
    # main000()

    current_densities, voltages = main()
    V_currentdensity(current_densities, voltages)

    # H2_Current()





