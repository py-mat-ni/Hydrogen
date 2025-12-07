from idaes.core import Component, LiquidPhase, PhaseType
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.state_definitions import FpcTP
from idaes.models.properties.modular_properties.pure import Perrys, RPP4, NIST
from pyomo.environ import units as pyunits, Param


# --- 自定义 KOH 参数类 ---
class ConstantKOHRho:
    @staticmethod
    def build_parameters(cobj):
        cobj.const_liq_rho = Param(
            initialize=22812,
            mutable=True,
            units=pyunits.mol / pyunits.m ** 3,
            doc="Constant molar density for KOH(aq)"
        )

    @staticmethod
    def return_expression(b, cobj, T):
        return cobj.const_liq_rho


class ConstantKOHEnthalpy:
    @staticmethod
    def build_parameters(cobj):
        cobj.const_liq_cp = Param(
            initialize=168.0,
            mutable=True,
            units=pyunits.J / pyunits.mol / pyunits.K,
            doc="Constant molar Cp for KOH(aq)"
        )
        cobj.enth_mol_form_liq_comp_ref = Param(
            initialize=-482.4e3,
            mutable=True,
            units=pyunits.J / pyunits.mol,
            doc="Standard formation enthalpy for KOH(liq)"
        )

    @staticmethod
    def return_expression(b, cobj, T):
        T_ref = b.params.temperature_ref
        H_ref = cobj.enth_mol_form_liq_comp_ref
        Cp = cobj.const_liq_cp
        return H_ref + Cp * (T - T_ref)

class NonVolatilePsat:
    @staticmethod
    def build_parameters(cobj):
        pass
    @staticmethod
    def return_expression(b, cobj, T):
        # 1e-25 Pa ≈ 10^-30 bar，绝对不可能挥发
        return 1e-25 * pyunits.Pa
# --- 液相配置字典 ---
configuration = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },

    "components": {
        # ============================================================
        #  H2O (Liquid)
        # ============================================================
        "H2O": {
            "type": Component,
            "elemental_composition": {"H": 2, "O": 1},
            "valid_phase_types": (PhaseType.liquidPhase,),  # 仅液相

            # 液相属性
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "dens_mol_liq_comp": Perrys,

            # ★★★ 关键：液相必须包含 Psat 计算，用于平衡计算 (Raoult定律)
            "pressure_sat_comp": NIST,

            "parameter_data": {
                "mw": (18.0153e-3, pyunits.kg / pyunits.mol),
                "pressure_crit": (220.64e5, pyunits.Pa),
                "temperature_crit": (647.14, pyunits.K),

                "cp_mol_liq_comp_coeff": {
                    '1': (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                    '2': (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K ** 2),
                    '3': (8.125, pyunits.J / pyunits.kmol / pyunits.K ** 3),
                    '4': (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K ** 4),
                    '5': (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K ** 5),
                },

                "enth_mol_form_liq_comp_ref": (-285.83e3, pyunits.J / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (69.95, pyunits.J / pyunits.mol / pyunits.K),

                "dens_mol_liq_comp_coeff": {
                    'eqn_type': 2,
                    '1': (4.9669, pyunits.kmol / pyunits.m ** 3),
                    '2': (2.7788e-1, pyunits.kmol / pyunits.m ** 3 / pyunits.K),
                    '3': (6.4713e2, pyunits.kmol / pyunits.m ** 3 / pyunits.K ** 2),
                    '4': (1.874e-1, pyunits.kmol / pyunits.m ** 3 / pyunits.K ** 3),
                },

                "pressure_sat_comp_coeff": {
                    'A': (4.6543, None),
                    'B': (1453.264, pyunits.K),
                    'C': (-64.848, pyunits.K),
                },
            },
        },

        # ============================================================
        #  KOH (Liquid)
        # ============================================================
        "KOH": {
            "type": Component,
            "valid_phase_types": (PhaseType.liquidPhase,),
            "elemental_composition": {"K": 1, "O": 1, "H": 1},

            "dens_mol_liq_comp": ConstantKOHRho,
            "enth_mol_liq_comp": ConstantKOHEnthalpy,
            "pressure_sat_comp": NonVolatilePsat, # KOH 不可能挥发
            "parameter_data": {
                "mw": (56.1056e-3, pyunits.kg / pyunits.mol),
            },
        },
    },

    "phases": {
        "Liq": {"type": LiquidPhase, "equation_of_state": Ideal},
    },
    "state_definition": FpcTP,
    "state_bounds": {
        "flow_mol_phase_comp": (1e-20, 1, 200, pyunits.mol / pyunits.s),
        "temperature": (1e-10, 353.15, 1000, pyunits.K),
        "pressure": (1e-10, 5e5, 1000e5, pyunits.Pa)
    },
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
}





