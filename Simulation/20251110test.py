import pyomo.environ as pyo
from pyomo.environ import ConcreteModel, SolverFactory, value, units as pyunits
from idaes.core import FlowsheetBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,)
# 导入你提供的 GasWasher 模型
from creating_model import GasWasher
from pyomo.environ import Var
from idaes.core.util import DiagnosticsToolbox
import pyomo.environ as pyo
from pyomo.environ import value, units as pyunits
from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import calculate_scaling_factors


import Liq as liq
import Vap as gas



def main():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # ==== 使用独立的物性包 ====
    m.fs.gas_props = GenericParameterBlock(**gas.configuration)
    m.fs.liq_props = GenericParameterBlock(**liq.configuration)

    # ==== 创建 GasWasher ====
    m.fs.washer = GasWasher(
        gas_property_package=m.fs.gas_props,
        liq_property_package=m.fs.liq_props,
    )

    # ==== 设置入口条件（略） ====
    t0 = m.fs.time.first()
    print("Initial DOF =", degrees_of_freedom(m))
    gm = m.fs.washer.gas_cv.properties_in[t0]
    print("Initial DOF =", degrees_of_freedom(m))
    gm.flow_mol_phase_comp["Vap", "H2"].fix(16)
    gm.flow_mol_phase_comp["Vap", "O2"].fix(4)
    gm.flow_mol_phase_comp["Vap", "H2O"].fix(1)
    print("Initial DOF =", degrees_of_freedom(m))
    gm.temperature.fix(350)
    gm.pressure.fix(1e5)
    print("Initial DOF =", degrees_of_freedom(m))

    lm = m.fs.washer.liq_cv.properties_in[t0]
    print("Initial DOF =", degrees_of_freedom(m))
    lm.flow_mol_phase_comp["Liq", "H2O"].fix(100)
    lm.flow_mol_phase_comp["Liq", "KOH"].fix(10)
    print("Initial DOF =", degrees_of_freedom(m))
    lm.temperature.fix(298.15)
    lm.pressure.fix(1e5)
    print("Initial DOF =", degrees_of_freedom(m))

    m.fs.washer.KGa_H2O.fix(10000)
    m.fs.washer.U_Area.fix(1000)
    m.fs.washer.gas_outlet.pressure.fix(1e5)
    print("Initial DOF =", degrees_of_freedom(m))

    # ==== 初始化 ====
    m.fs.washer.initialize()

    # ==== 求解 ====
    solver = SolverFactory("ipopt")
    solver.solve(m, tee=True)

    # 8) 打印一部分结果
    print("\n===== Results =====")
    print("Gas outlet T  = ", value(m.fs.washer.gas_cv.properties_out[t0].temperature), "K")
    print("Liq outlet T  = ", value(m.fs.washer.liq_cv.properties_out[t0].temperature), "K")
    print("Gas outlet y_H2O = ",
          value(m.fs.washer.gas_cv.properties_out[t0].mole_frac_phase_comp["Vap", "H2O"]))
    print("Liq outlet x_H2O = ",
          value(m.fs.washer.liq_cv.properties_out[t0].mole_frac_phase_comp["Liq", "H2O"]))
    print("N_H2O (mol/s) = ", value(m.fs.washer.N_H2O[t0]))
    print("Q_heat (W)    = ", value(m.fs.washer.Q_heat[t0]))
    import os, time
    if os.path.exists("Hydro.json"):
        os.remove("Hydro.json")
    m.fs.visualize("Hydro")
    while True:
        time.sleep(1)
    # from idaes.core.util import DiagnosticsToolbox
    # dt = DiagnosticsToolbox(m)
    # dt.report_structural_issues()
    # dt.display_components_with_inconsistent_units()
    # dt.display_underconstrained_set()
    # dt.display_overconstrained_set()
    # dt.display_potential_evaluation_errors()

if __name__ == "__main__":
    main()







