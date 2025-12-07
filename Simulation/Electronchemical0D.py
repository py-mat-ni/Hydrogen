from pyomo.environ import (
    Var, Param, Constraint, units as pyunits, exp, log, value
)
from idaes.core import (
    UnitModelBlockData, declare_process_block_class, useDefault
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.constants import Constants
from idaes.models.unit_models import Separator
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    ControlVolume0DBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    useDefault,
)
from idaes.core.util import scaling as iscale
from pyomo.common.config import ConfigValue, ConfigBlock, In, Bool
from pyomo.environ import (
    Var,
    Constraint,
    Expression,
    NonNegativeReals,
    units as pyunits,
)
# Electronchemical0D.py - 完美修正版 (含压力平衡)
from pyomo.environ import (
    Var, Param, Constraint, units as pyunits, exp, log, value
)
import pyomo as pyo
from idaes.core import (
    UnitModelBlockData, declare_process_block_class, useDefault,
    EnergyBalanceType
)
from idaes.core import ControlVolume0DBlock
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block
)
from pyomo.common.config import ConfigValue, ConfigDict
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor
from idaes.core.util.scaling import constraint_scaling_transform
from idaes.core.util.exceptions import InitializationError
import idaes.logger as idaeslog
import pyomo.environ as pyo
from pyomo.core.base.reference import Reference

@declare_process_block_class("AlkalineElectrolyzer0D")
class AlkalineElectrolyzer0DData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume"
    ))

    CONFIG.declare("reaction_package", ConfigValue(
        default=None,
        domain=is_reaction_parameter_block,
        description="Reaction package for electrolysis"
    ))
    CONFIG.declare("has_heat_transfer", ConfigValue(
        default=False,  # 默认为 False (绝热模式)
        domain=Bool,
        description="Set to True to activate heat transfer term"
    ))


    def build(self):
        super().build()

        self.control_volume = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.property_package,
            reaction_package=self.config.reaction_package
        )

        # 2. 添加状态块
        self.control_volume.add_state_blocks(has_phase_equilibrium=False)

        # 3. 添加反应块
        self.control_volume.add_reaction_blocks(has_equilibrium=False)

        # 4. 添加能量平衡
        self.control_volume.add_energy_balances(
            balance_type=EnergyBalanceType.enthalpyTotal,
            has_heat_transfer=self.config.has_heat_transfer,
            has_work_transfer=True
        )
        self.control_volume.add_material_balances(
            balance_type=MaterialBalanceType.componentTotal,
            has_rate_reactions=True
        )


        # 5. 添加端口
        self.add_inlet_port()
        self.add_outlet_port()
        t_set = self.flowsheet().time
        if self.config.has_heat_transfer:
            self.heat_duty = Reference(self.control_volume.heat)
        # 停用 ReactionBlock 内部的 Arrhenius 计算
        for t in t_set:
            if hasattr(self.control_volume.reactions[t], "reaction_rate_eq"):
                self.control_volume.reactions[t].reaction_rate_eq.deactivate()
            if hasattr(self.control_volume.reactions[t], "equilibrium_constraint"):
                self.control_volume.reactions[t].equilibrium_constraint.deactivate()


        self.current = Var(t_set, initialize=5000, bounds=(1e-6, 20000),units=pyunits.ampere, doc="Total Current")
        self.voltage = Var(t_set, initialize=1.8, bounds=(1.0, 3.5),units=pyunits.volt, doc="Cell Voltage")
        self.power = Var(t_set, initialize=9000, bounds=(1e-6, None),units=pyunits.watt, doc="Total Power")
        self.current_density = Var(t_set, initialize=4000, bounds=(1e-6, 10000),units=pyunits.ampere / pyunits.m ** 2)

        # 参数
        self.number_cells = Param(initialize=326, mutable=True)
        self.area_cell = Param(initialize=2.66, units=pyunits.m ** 2)
        self.faraday_efficiency = Param(initialize=0.86, mutable=True)
        self.F = Param(initialize=96485.33, units=pyunits.coulomb / pyunits.mol)
        self.z = Param(initialize=2)

        # U-I 参数
        self.alpha_1 = Param(initialize=0.8, units=pyunits.ohm * pyunits.cm ** 2)
        self.alpha_2 = Param(initialize=-0.00763, units=pyunits.ohm * pyunits.cm ** 2 / pyunits.K)
        self.s_coeff = Param(initialize=0.1795, units=pyunits.volt)
        self.beta_1 = Param(initialize=20.0, units=pyunits.cm ** 2 / pyunits.ampere)
        self.beta_2 = Param(initialize=0.1, units=pyunits.cm ** 2 * pyunits.K / pyunits.ampere)
        self.beta_3 = Param(initialize=3.5e5, units=pyunits.cm ** 2 * pyunits.K ** 2 / pyunits.ampere)


        @self.Constraint(t_set)
        def eq_faraday_driver(b, t):
            rxn_extent = b.control_volume.rate_reaction_extent[t, "water_splitting"] / (pyunits.mol / pyunits.s)
            I = b.current[t] / pyunits.ampere
            N = b.number_cells
            eff = b.faraday_efficiency
            z = b.z
            F = b.F / (pyunits.coulomb / pyunits.mol)

            prod_rate = (N * I * eff) / (z * F)
            return rxn_extent == prod_rate

        self.Pv_H2O_bar = Var(t_set, initialize=0.47, bounds=(1e-6, 100),units=pyunits.bar)
        self.Pv_KOH = Var(t_set, initialize=35000, bounds=(1e-6, 1e7),units=pyunits.Pa)
        self.activity_water = Var(t_set, initialize=0.7, bounds=(1e-6, 1.0),units=pyunits.dimensionless)
        self.coef_a = Var(initialize=-0.1, units=pyunits.dimensionless)
        self.coef_b = Var(initialize=0.9, units=pyunits.dimensionless)
        self.molality = Param(initialize=5.9407, units=pyunits.mol / pyunits.kg)

        @self.Constraint(t_set)
        def eq_Pv_H2O_pure(b, t):
            T_degC = (b.control_volume.properties_out[t].temperature / pyunits.K) - 273.15
            P_bar = b.Pv_H2O_bar[t] / pyunits.bar
            return log(P_bar) / log(10) == 5.1962 - 1730.63 / (T_degC + 233.426)

        @self.Constraint()
        def eq_coef_a(b):
            M = b.molality / (pyunits.mol / pyunits.kg)  # 提取无量纲数值
            return b.coef_a == -0.0151 * M - 1.6788e-3 * (M ** 2) + 2.2588e-5 * (M ** 3)

        @self.Constraint()
        def eq_coef_b(b):
            M = b.molality / (pyunits.mol / pyunits.kg)
            return b.coef_b == 1 - 1.2062e-3 * M + 5.6024e-4 * (M ** 2) - 7.8228e-6 * (M ** 3)

        @self.Constraint(t_set)
        def eq_Pv_KOH_solution(b, t):
            P_pure_bar = b.Pv_H2O_bar[t] / pyunits.bar
            P_koh_pa = b.Pv_KOH[t] / pyunits.Pa
            term = 2.302 * b.coef_a + b.coef_b * log(P_pure_bar)
            return P_koh_pa == exp(term) * 1e5

        @self.Constraint(t_set)
        def eq_water_activity(b, t):
            M = b.molality / (pyunits.mol / pyunits.kg)
            T_val = b.control_volume.properties_out[t].temperature / pyunits.K

            exponent = -0.05192 * M + 0.003302 * (M ** 2) + (3.177 * M - 2.131 * (M ** 2)) / T_val

            return b.activity_water[t] == exp(exponent)



        # 假设无压降: P_out = P_in
        @self.Constraint(t_set)
        def eq_pressure_balance(b, t):
            return b.control_volume.properties_out[t].pressure == b.control_volume.properties_in[t].pressure



        # (D) 电化学约束
        @self.Constraint(t_set)
        def eq_current_density(b, t):
            return b.current_density[t] == b.current[t] / b.area_cell

        @self.Constraint(t_set)
        def eq_voltage(b, t):
            # 温度 (K)
            T_val = b.control_volume.properties_out[t].temperature / pyunits.K

            # 压力 (Pa)
            P_sys_val = b.control_volume.properties_out[t].pressure / pyunits.Pa
            P_koh_val = b.Pv_KOH[t] / pyunits.Pa


            i_val = (b.current_density[t] / (pyunits.ampere / pyunits.m ** 2)) * 1e-4

            a_water_val = b.activity_water[t]


            alpha1 = b.alpha_1 / (pyunits.ohm * pyunits.cm ** 2)
            alpha2 = b.alpha_2 / (pyunits.ohm * pyunits.cm ** 2 / pyunits.K)
            s = b.s_coeff / pyunits.volt
            beta1 = b.beta_1 / (pyunits.cm ** 2 / pyunits.ampere)
            beta2 = b.beta_2 / (pyunits.cm ** 2 * pyunits.K / pyunits.ampere)
            beta3 = b.beta_3 / (pyunits.cm ** 2 * pyunits.K ** 2 / pyunits.ampere)

            R = 8.31446
            F = 96485.33
            z = 2.0
            U0 = 1.5184 - 1.5421e-3 * T_val + 9.526e-5 * T_val * log(T_val) + 9.84e-8 * T_val ** 2
            P_ratio = (P_sys_val - P_koh_val) / 1e5
            term_nernst = ((R * T_val) / (z * F)) * log((P_ratio ** 1.5) / a_water_val)
            U_rev = U0 + term_nernst
            T_C = T_val - 273.15
            r_ohm = alpha1 + alpha2 * T_C
            U_ohm = r_ohm * i_val
            beta = beta1 + beta2 / T_C + beta3 / (T_C ** 2)
            arg_log = beta * i_val + 1
            U_act = s * log(arg_log)
            U_total = U_rev + U_ohm + U_act
            return b.voltage[t] / pyunits.volt == U_total

        @self.Constraint(t_set)
        def eq_power(b, t):
            return b.power[t] == b.voltage[t] * b.current[t] * b.number_cells

        @self.Constraint(t_set)
        def eq_energy_work(b, t):
            return b.control_volume.work[t] == b.power[t]

        @self.Constraint(t_set)
        def eq_water_evaporation_balance(b, t):
            P_vap = b.Pv_KOH[t] / pyunits.Pa
            P_sys = b.control_volume.properties_out[t].pressure / pyunits.Pa

            flow_H2 = b.control_volume.properties_out[t].flow_mol_phase_comp["Vap", "H2"]
            flow_O2 = b.control_volume.properties_out[t].flow_mol_phase_comp["Vap", "O2"]
            flow_H2O = b.control_volume.properties_out[t].flow_mol_phase_comp["Vap", "H2O"]

            # 使用无量纲数值乘法形式
            return flow_H2O * (P_sys - P_vap) == (flow_H2 + flow_O2) * P_vap




