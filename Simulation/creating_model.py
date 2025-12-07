# gas_washer_example.py

from pyomo.environ import (
    ConcreteModel, SolverFactory, value, units as pyunits,
    Var, NonNegativeReals, Constraint, Expression
)
from pyomo.common.config import ConfigBlock, ConfigValue, In
from idaes.core import (
    FlowsheetBlock, declare_process_block_class,
    UnitModelBlockData, ControlVolume0DBlock,
    MaterialBalanceType, EnergyBalanceType, MomentumBalanceType,
    PhysicalParameterBlock, Component, VaporPhase, LiquidPhase, PhaseType
)
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block
)
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock
)
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.state_definitions import FpcTP
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import IdealBubbleDew
from idaes.models.properties.modular_properties.pure import RPP4, NIST
from idaes.models.properties.modular_properties.pure.Perrys import Perrys
from pyomo.environ import Param

# =============================
#   GasWasher 单元模型定义
# =============================
@declare_process_block_class("GasWasher")
class GasWasherData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "gas_property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            doc="Gas phase property package",
        ),
    )
    CONFIG.declare(
        "gas_property_package_args",
        ConfigBlock(implicit=True, description="Arguments for gas props"),
    )
    CONFIG.declare(
        "liq_property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            doc="Liquid phase property package",
        ),
    )
    CONFIG.declare(
        "liq_property_package_args",
        ConfigBlock(implicit=True, description="Arguments for liquid props"),
    )

    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.componentPhase,
            domain=In(MaterialBalanceType),
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.enthalpyTotal,
            domain=In(EnergyBalanceType),
        ),
    )
    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.none,
            domain=In(MomentumBalanceType),
        ),
    )

    def build(self):
        super().build()
        tset = self.flowsheet().time

        # --- 气相控制体 ---
        self.gas_cv = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=False,
            property_package=self.config.gas_property_package,
            property_package_args=self.config.gas_property_package_args,
        )
        self.gas_cv.add_state_blocks(has_phase_equilibrium=False)
        self.gas_cv.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=True,
        )
        self.gas_cv.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_enthalpy_transfer=True,
        )

        # --- 液相控制体 ---
        self.liq_cv = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=False,
            property_package=self.config.liq_property_package,
            property_package_args=self.config.liq_property_package_args,
        )
        self.liq_cv.add_state_blocks(has_phase_equilibrium=False)
        self.liq_cv.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=True,
        )
        self.liq_cv.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_enthalpy_transfer=True,
        )

        # --- 端口 ---
        self.add_inlet_port("gas_inlet", self.gas_cv)
        self.add_outlet_port("gas_outlet", self.gas_cv)
        self.add_inlet_port("liquid_inlet", self.liq_cv)
        self.add_outlet_port("liquid_outlet", self.liq_cv)

        # --- 传质/传热参数 ---
        self.KGa_H2O = Var(
            initialize=0.1,
            domain=NonNegativeReals,
            units=pyunits.mol / pyunits.s,
            doc="Overall mass transfer coefficient * area for H2O",
        )
        self.U_Area = Var(
            initialize=1000.0,
            domain=NonNegativeReals,
            units=pyunits.W / pyunits.K,
            doc="Overall heat transfer coefficient * area",
        )

        self.N_H2O = Var(tset, units=pyunits.mol / pyunits.s, doc="H2O transfer rate")
        self.Q_heat = Var(tset, units=pyunits.W, doc="Heat transfer")

        # --- 传质驱动力 ---
        @self.Expression(tset)
        def y_H2O_bulk(b, t):
            return b.gas_cv.properties_out[t].mole_frac_phase_comp["Vap", "H2O"]

        @self.Expression(tset)
        def y_H2O_star(b, t):
            # Raoult 定律： y* = x_H2O * Psat / P_gas
            Psat = b.liq_cv.properties_out[t].pressure_sat_comp["H2O"]
            x_H2O = b.liq_cv.properties_out[t].mole_frac_phase_comp["Liq", "H2O"]
            P_g = b.gas_cv.properties_out[t].pressure
            return (x_H2O * Psat) / P_g

        @self.Expression(tset)
        def driving_force_H2O(b, t):
            return b.y_H2O_bulk[t] - b.y_H2O_star[t]

        # --- 传质速率方程 ---
        @self.Constraint(tset)
        def eq_MT_rate(b, t):
            return b.N_H2O[t] == b.KGa_H2O * b.driving_force_H2O[t]

        # --- 传热方程 ---
        @self.Constraint(tset)
        def eq_HT_rate(b, t):
            T_g = b.gas_cv.properties_out[t].temperature
            T_l = b.liq_cv.properties_out[t].temperature
            return b.Q_heat[t] == b.U_Area * (T_g - T_l)

        # --- 质量传递映射（H2O） ---
        @self.Constraint(tset)
        def eq_gas_MT_H2O(b, t):
            return b.gas_cv.mass_transfer_term[t, "Vap", "H2O"] == -b.N_H2O[t]

        @self.Constraint(tset)
        def eq_liq_MT_H2O(b, t):
            return b.liq_cv.mass_transfer_term[t, "Liq", "H2O"] == b.N_H2O[t]

        # --- 气相惰性组分（H2, O2）不发生传质 ---
        gas_comps = self.config.gas_property_package.component_list

        @self.Constraint(tset, gas_comps)
        def eq_gas_MT_inerts(b, t, j):
            if j == "H2O":
                return Constraint.Skip
            return b.gas_cv.mass_transfer_term[t, "Vap", j] == 0

        # --- 液相惰性组分（KOH）不发生传质 ---
        liq_comps = self.config.liq_property_package.component_list

        @self.Constraint(tset, liq_comps)
        def eq_liq_MT_inerts(b, t, j):
            if j == "H2O":
                return Constraint.Skip
            return b.liq_cv.mass_transfer_term[t, "Liq", j] == 0

        # --- 焓传递（相变 + 显热） ---
        @self.Expression(tset)
        def enth_transfer_flow(b, t):
            h_vap = b.gas_cv.properties_out[t].enth_mol_phase_comp["Vap", "H2O"]
            return b.N_H2O[t] * h_vap

        @self.Constraint(tset)
        def eq_gas_enthalpy_balance(b, t):
            return b.gas_cv.enthalpy_transfer[t] == -b.enth_transfer_flow[t] - b.Q_heat[t]

        @self.Constraint(tset)
        def eq_liq_enthalpy_balance(b, t):
            return b.liq_cv.enthalpy_transfer[t] == b.enth_transfer_flow[t] + b.Q_heat[t]

        # --- 气液同压 ---
        @self.Constraint(tset)
        def eq_P(b, t):
            return b.gas_cv.properties_out[t].pressure == b.liq_cv.properties_out[t].pressure

    # ====== 自定义初始化 ======
    def initialize_build(self, *args, **kwargs):
        # 简单做法：先把 N_H2O, Q_heat 固定为 0，初始化两个 CV，再解开
        t0 = self.flowsheet().time.first()

        self.N_H2O[t0].fix(0.0)
        self.Q_heat[t0].fix(0.0)

        self.gas_cv.initialize(*args, **kwargs)
        self.liq_cv.initialize(*args, **kwargs)

        self.N_H2O[t0].unfix()
        self.Q_heat[t0].unfix()

        return {}
