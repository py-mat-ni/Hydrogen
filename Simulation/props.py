from idaes.core import Component, LiquidPhase, VaporPhase, PhaseType
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import IdealBubbleDew
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity
from idaes.core import PhysicalParameterBlock
from idaes.core.util.misc import add_object_reference
from idaes.core import declare_process_block_class
from pyomo.environ import units as pyunits
from idaes.models.properties.modular_properties.pure.RPP4 import RPP4

# --- 纯物性相关式（内置可用的） ---
from idaes.models.properties.modular_properties.pure.RPP4 import RPP4  # ideal-gas Cp, Psat forms
from idaes.models.properties.modular_properties.pure.Perrys import Perrys  # liquid Cp, rho

# -------【可选】给“KOH 液相专用：常值密度/热容”的轻量实现 -----------
# 如果你没有直接的KOH(aq)密度/比热经验式，可用这个“常值版本”先跑通。
# 方法：借用 Perrys 接口风格的“自定义类”写出恒定返回值（简单而粗糙）。

from pyomo.environ import Param, value

# === props.py 中 ===

def ConstantLiquidCp_build_parameters(cobj):
    # 定义常值 Cp 参数（J/kmol/K）
    cobj.const_liq_cp = Param(
        initialize=80e3,
        mutable=True,
        units=pyunits.J/pyunits.kmol/pyunits.K,
        doc="Constant liquid molar Cp for KOH(aq)"
    )


def ConstantLiquidCp_func(b, cobj, T):
    # 返回常值 Cp，不依赖温度
    return value(cobj.const_liq_cp)


def ConstantLiquidRho_build_parameters(cobj):
    # 定义常值密度参数（kmol/m3）
    cobj.const_liq_rho = Param(
        initialize=60.0,
        mutable=True,
        units=pyunits.kmol/pyunits.m**3,
        doc="Constant liquid molar density for KOH(aq)"
    )


def ConstantLiquidRho_func(b, cobj, T):
    return value(cobj.const_liq_rho)

# ----------------------------------------------------------

configuration = {
    # 1) 组分
    "components": {
        # 氢气（只在气相参与；若你日后考虑溶解，可加 Henry 常数与液相有效性）
        "H2": {
            "type": Component,
            "valid_phase_types": (PhaseType.vaporPhase,),  # 仅气相
            "elemental_composition": {"H": 2},
            "enth_mol_ig_comp": RPP4,            # 理想气热容/焓
            # 氢气无 Psat / 无液相纯物性
            "parameter_data": {
                "mw": (2.01588e-3, pyunits.kg/pyunits.mol),
                "cp_mol_ig_comp_coeff": {
                    "A": (29.11, pyunits.J / pyunits.mol / pyunits.K),
                    "B": ((-0.1916e-2), pyunits.J / pyunits.mol / pyunits.K ** 2),
                    "C": ((0.4003e-5), pyunits.J / pyunits.mol / pyunits.K ** 3),
                    "D": ((-0.8704e-9), pyunits.J / pyunits.mol / pyunits.K ** 4),
                },
                "pressure_crit": (12.964e5, pyunits.Pa),
                "temperature_crit": (33.19, pyunits.K),
                "enth_mol_form_vap_comp_ref": (0, pyunits.J/pyunits.mol),
            },
        },
        # 氧气（同氢气）
        "O2": {
            "type": Component,
            "valid_phase_types": (PhaseType.vaporPhase,),
            "elemental_composition": {"O": 2},
            "enth_mol_ig_comp": RPP4,
            "parameter_data": {
                "mw": (31.998e-3, pyunits.kg/pyunits.mol),
                "pressure_crit": (50.43e5, pyunits.Pa),
                "temperature_crit": (154.58, pyunits.K),
                "enth_mol_form_vap_comp_ref": (0, pyunits.J/pyunits.mol),
                "cp_mol_ig_comp_coeff": {
                    "A": (25.48, pyunits.J / pyunits.mol / pyunits.K),
                    "B": ((1.520e-2), pyunits.J / pyunits.mol / pyunits.K ** 2),
                    "C": ((-0.7155e-5), pyunits.J / pyunits.mol / pyunits.K ** 3),
                    "D": ((1.312e-9), pyunits.J / pyunits.mol / pyunits.K ** 4),
                },
            },
        },
        # 水（既有气相也有液相；参与 VLE）
        "H2O": {
            "type": Component,
            "elemental_composition": {"H": 2, "O": 1},
            "enth_mol_ig_comp": RPP4,            # 气相 Cp/焓（RPP4）
            "pressure_sat_comp": RPP4,           # Psat（RPP4）
            "enth_mol_liq_comp": Perrys,         # 液相 Cp/焓（Perrys）
            "dens_mol_liq_comp": Perrys,         # 液相密度（Perrys）
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},  # 对水做 VLE 逸度平衡
            "parameter_data": {
                "mw": (18.01528e-3, pyunits.kg/pyunits.mol),
                "pressure_crit": (22.064e6, pyunits.Pa),
                "temperature_crit": (647.14, pyunits.K),
                "enth_mol_form_vap_comp_ref": (-241.826e3, pyunits.J/pyunits.mol),
                "enth_mol_form_liq_comp_ref": (-285.83e3, pyunits.J/pyunits.mol),

                "cp_mol_ig_comp_coeff": {
                    "A": (32.24, pyunits.J / pyunits.mol / pyunits.K),
                    "B": (0.001924, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (1.055e-5, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (-3.596e-9, pyunits.J / pyunits.mol / pyunits.K**4),
                },
        # ---- 液相Cp (Perry’s 5项式) ----
                "cp_mol_liq_comp_coeff": {
                    "1": (2.7637e2, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2.0901e-1, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (8.125e-4, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (-1.4116e-6, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (9.3701e-10, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
        # ---- 液相密度 (Perry’s Type 1) ----
                "dens_mol_liq_comp_coeff": {
                    "eqn_type": 1,
                    "1": (1.0439, pyunits.kmol / pyunits.m**3),
                    "2": (0.2362, None),
                    "3": (647.13, pyunits.K),
                    "4": (0.6366, None),
                },
        # ---- 饱和蒸汽压 (Antoine/RPP4型) ----
                "pressure_sat_comp_coeff": {
                    "A": (-7.85951783, None),
                    "B": (1.84408259, None),
                    "C": (-11.7866497, None),
                    "D": (22.6807411, None),
                },
            },
        },


        # KOH：仅液相、不挥发（不参与 VLE）
        "KOH": {
            "type": Component,
            "valid_phase_types": (PhaseType.liquidPhase,),  # 仅液相
            "elemental_composition": {"K": 1, "O": 1, "H": 1},
            # 自定义：用常值密度 / 常值Cp（你也可以换成自己的经验式类）
            "dens_mol_liq_comp": ConstantLiquidRho_func,
            "enth_mol_liq_comp": ConstantLiquidCp_func,
            # 不提供 Psat、不提供气相 Cp（表示不考虑其挥发/蒸汽相）
            "parameter_data": {
                "mw": (56.1056e-3, pyunits.kg/pyunits.mol),
                "enth_mol_form_vap_comp_ref": (0, pyunits.J/pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "eqn_type": 1,
                    "1": (61.0, pyunits.kmol / pyunits.m ** 3),
                    "2": (0.0, None),
                    "3": (298.15, pyunits.K),
                    "4": (0.0, None),
                },
                "cp_mol_liq_comp_coeff": {
                    "1": (100e3, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (0, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (0, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
            },
        },
    },


    # 2) 相
    "phases": {
        "Vap": {"type": VaporPhase, "equation_of_state": Ideal},   # 理想气
        "Liq": {"type": LiquidPhase, "equation_of_state": Ideal},  # 理想液（做理想混合）
    },

    # 3) 基本单位
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },

    # 4) 状态定义与边界
    "state_definition": FTPx,
    "state_bounds": {
        "flow_mol": (0, 10, 1e3, pyunits.mol/pyunits.s),
        "temperature": (273.15, 350, 600, pyunits.K),
        "pressure": (5e4, 1e5, 5e6, pyunits.Pa),
    },
    "pressure_ref": (1e5, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),

    # 5) 相平衡（仅对水）
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": IdealBubbleDew,

    # 6) 缩放示例（按需加）
    "default_scaling_factors": {
        "mole_frac_comp[hydrogen]": 1,
        "mole_frac_comp[oxygen]": 1,
        "mole_frac_comp[water]": 1,
        "mole_frac_comp[KOH]": 1,
    },

}






