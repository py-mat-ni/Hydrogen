from idaes.core import Component, LiquidPhase, VaporPhase, PhaseType
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.state_definitions import FpcTP
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import IdealBubbleDew
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity
from idaes.models.properties.modular_properties.pure import RPP4, NIST, Perrys
from pyomo.environ import units as pyunits, Param


# --- 自定义常值 Cp / 密度简化实现 ---
class ConstantKOHRho:
    @staticmethod
    def build_parameters(cobj):
        cobj.const_liq_rho = Param(
            initialize=22812,  # 22.812 kmol/m3 = 22812 mol/m3
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
            initialize=168.0,  # J/mol/K
            mutable=True,
            units=pyunits.J / pyunits.mol / pyunits.K,
            doc="Constant molar Cp for KOH(aq)"
        )

        cobj.enth_mol_form_liq_comp_ref = Param(
            initialize=-482.4e3,  # J/mol
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


# --- 配置字典 ---
configuration_ELY = {

    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },

    "components": {
        # ============================================================
        #  WATER
        # ============================================================
        "H2O": {
            "type": Component,
            "elemental_composition": {"H": 2, "O": 1},
            "valid_phase_types": (PhaseType.vaporPhase, PhaseType.liquidPhase),

            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            "entr_mol_ig_comp": NIST,

            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "dens_mol_liq_comp": Perrys,

            "pressure_sat_comp": NIST,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},

            "parameter_data": {
                "mw": (18.0153e-3, pyunits.kg / pyunits.mol),

                "pressure_crit": (220.64e5, pyunits.Pa),
                "temperature_crit": (647.14, pyunits.K),

                "cp_mol_ig_comp_coeff": {
                    'A': (30.092, pyunits.J / (pyunits.mol * pyunits.K)),
                    'B': (6.832514, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK)),
                    'C': (6.793435, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK ** 2)),
                    'D': (-2.53448, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK ** 3)),
                    'E': (0.082139, pyunits.J * pyunits.kiloK ** 2 / (pyunits.mol * pyunits.K)),
                    'F': (-250.881, pyunits.kiloJ / pyunits.mol),
                    'G': (223.3967, pyunits.J / (pyunits.mol * pyunits.K)),
                    'H': (-241.8264, pyunits.kiloJ / pyunits.mol),
                },

                "cp_mol_liq_comp_coeff": {
                    '1': (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                    '2': (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K ** 2),
                    '3': (8.125, pyunits.J / pyunits.kmol / pyunits.K ** 3),
                    '4': (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K ** 4),
                    '5': (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K ** 5),
                },

                "enth_mol_form_liq_comp_ref": (-285.83e3, pyunits.J / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (69.95, pyunits.J / pyunits.mol / pyunits.K),

                "pressure_sat_comp_coeff": {
                    'A': (4.6543, None),
                    'B': (1453.264, pyunits.K),
                    'C': (-64.848, pyunits.K),
                },

                "dens_mol_liq_comp_coeff": {
                    'eqn_type': 2,
                    '1': (4.9669, pyunits.kmol / pyunits.m ** 3),
                    '2': (2.7788e-1, pyunits.kmol / pyunits.m ** 3 / pyunits.K),
                    '3': (6.4713e2, pyunits.kmol / pyunits.m ** 3 / pyunits.K ** 2),
                    '4': (1.874e-1, pyunits.kmol / pyunits.m ** 3 / pyunits.K ** 3),
                },
            },
        },

        # ============================================================
        #  H2 (gas only) - 修正 pressure_sat_comp
        # ============================================================
        "H2": {
            "type": Component,
            "elemental_composition": {"H": 2},
            "valid_phase_types": (PhaseType.vaporPhase,),

            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            "entr_mol_ig_comp": NIST,

            "pressure_sat_comp": None,  # 修正：永久气体不计算饱和压
            "parameter_data": {
                "mw": (2.01588e-3, pyunits.kg / pyunits.mol),
                "pressure_crit": (13e5, pyunits.Pa),
                "temperature_crit": (33.18, pyunits.K),
                "cp_mol_ig_comp_coeff": {
                    'A': (33.066178, pyunits.J / (pyunits.mol * pyunits.K)),
                    'B': (-11.363417, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK)),
                    'C': (11.432816, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK ** 2)),
                    'D': (-2.772874, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK ** 3)),
                    'E': (-0.158558, (pyunits.J * pyunits.kiloK ** 2) / (pyunits.mol * pyunits.K)),
                    'F': (-9.980797, (pyunits.kiloJ / pyunits.mol)),
                    'G': (172.707974, pyunits.J / (pyunits.mol * pyunits.K)),
                    'H': (1e-10, (pyunits.kiloJ / pyunits.mol))
                },
            }
        },

        "O2": {
            "type": Component,
            "elemental_composition": {"O": 2},
            "valid_phase_types": (PhaseType.vaporPhase,),

            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            "entr_mol_ig_comp": NIST,

            "pressure_sat_comp": None,  # 永久气体不计算饱和压
            "parameter_data": {
                "mw": (31.9988e-3, pyunits.kg / pyunits.mol),
                "pressure_crit": (50.43e5, pyunits.Pa),
                "temperature_crit": (154.58, pyunits.K),
                "cp_mol_ig_comp_coeff": {
                    'A': (31.32234, pyunits.J / (pyunits.mol * pyunits.K)),
                    'B': (-20.23531, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK)),
                    'C': (57.86644, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK ** 2)),
                    'D': (-36.50624, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK ** 3)),
                    'E': (-0.007374, (pyunits.J * pyunits.kiloK ** 2) / (pyunits.mol * pyunits.K)),
                    'F': (-8.903471, (pyunits.kiloJ / pyunits.mol)),
                    'G': (246.7945, pyunits.J / (pyunits.mol * pyunits.K)),
                    'H': (1e-10, (pyunits.kiloJ / pyunits.mol))
                },
            }
        },

        # ============================================================
        #  KOH (liquid only)
        # ============================================================
        "KOH": {
            "type": Component,
            "valid_phase_types": (PhaseType.liquidPhase,),
            "elemental_composition": {"K": 1, "O": 1, "H": 1},

            "dens_mol_liq_comp": ConstantKOHRho,
            "enth_mol_liq_comp": ConstantKOHEnthalpy,

            "parameter_data": {
                "mw": (56.1056e-3, pyunits.kg / pyunits.mol),
            },
        },
    },

    "phases": {
        "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
        "Liq": {"type": LiquidPhase, "equation_of_state": Ideal},
    },

    # --- 相平衡全局配置 ---
    # "phases_in_equilibrium": [("Vap", "Liq")],
    # "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    # "bubble_dew_method": IdealBubbleDew,

    "state_definition": FpcTP,


    "state_bounds": {
        "flow_mol_phase_comp": (1e-20, 1, 20000, pyunits.mol / pyunits.s),
        "temperature": (300, 353.15, 400, pyunits.K),
        "pressure": (1e-10, 5e5, 1e7, pyunits.Pa)
    },
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
    # "default_scaling_factors": {
    #     ("temperature", None): 1e-2,
    #     ("pressure", None): 1e-5,     # 将 Pa 缩放至 bar 量级
    #
    #     # 流量缩放 (基于预期值)
    #     ("flow_mol_phase_comp", ("Liq", "H2O")): 1e-3, # ~800 mol/s
    #     ("flow_mol_phase_comp", ("Liq", "KOH")): 1e-2,
    #     ("flow_mol_phase_comp", ("Vap", "H2")): 1e-1,  # ~13 mol/s
    #     ("flow_mol_phase_comp", ("Vap", "O2")): 1e-1,
    #     ("flow_mol_phase_comp", ("Vap", "H2O")): 1,
    #
    #     # 关键：焓缩放。目标：(流量缩放) * (焓缩放) = (功率缩放 1e-6)
    #     # 1e-3 * 1e-3 = 1e-6
    #     ("enth_mol", None): 1e-3,
    #     ("enth_mol_phase", "Liq"): 1e-3,
    #     ("enth_mol_phase", "Vap"): 1e-3,
    # }
}











































configuration_VLE = {

    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },

    "components": {
        # ============================================================
        #  WATER
        # ============================================================
        "H2O": {
            "type": Component,
            "elemental_composition": {"H": 2, "O": 1},
            "valid_phase_types": (PhaseType.vaporPhase, PhaseType.liquidPhase),

            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            "entr_mol_ig_comp": NIST,

            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "dens_mol_liq_comp": Perrys,

            "pressure_sat_comp": NIST,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},

            "parameter_data": {
                "mw": (18.0153e-3, pyunits.kg / pyunits.mol),

                "pressure_crit": (220.64e5, pyunits.Pa),
                "temperature_crit": (647.14, pyunits.K),

                "cp_mol_ig_comp_coeff": {
                    'A': (30.092, pyunits.J / (pyunits.mol * pyunits.K)),
                    'B': (6.832514, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK)),
                    'C': (6.793435, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK ** 2)),
                    'D': (-2.53448, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK ** 3)),
                    'E': (0.082139, pyunits.J * pyunits.kiloK ** 2 / (pyunits.mol * pyunits.K)),
                    'F': (-250.881, pyunits.kiloJ / pyunits.mol),
                    'G': (223.3967, pyunits.J / (pyunits.mol * pyunits.K)),
                    'H': (-241.8264, pyunits.kiloJ / pyunits.mol),
                },

                "cp_mol_liq_comp_coeff": {
                    '1': (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                    '2': (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K ** 2),
                    '3': (8.125, pyunits.J / pyunits.kmol / pyunits.K ** 3),
                    '4': (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K ** 4),
                    '5': (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K ** 5),
                },

                "enth_mol_form_liq_comp_ref": (-285.83e3, pyunits.J / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (69.95, pyunits.J / pyunits.mol / pyunits.K),

                "pressure_sat_comp_coeff": {
                    'A': (4.6543, None),
                    'B': (1453.264, pyunits.K),
                    'C': (-64.848, pyunits.K),
                },

                "dens_mol_liq_comp_coeff": {
                    'eqn_type': 2,
                    '1': (4.9669, pyunits.kmol / pyunits.m ** 3),
                    '2': (2.7788e-1, pyunits.kmol / pyunits.m ** 3 / pyunits.K),
                    '3': (6.4713e2, pyunits.kmol / pyunits.m ** 3 / pyunits.K ** 2),
                    '4': (1.874e-1, pyunits.kmol / pyunits.m ** 3 / pyunits.K ** 3),
                },
            },
        },

        # ============================================================
        #  H2 (gas only) - 修正 pressure_sat_comp
        # ============================================================
        "H2": {
            "type": Component,
            "elemental_composition": {"H": 2},
            "valid_phase_types": (PhaseType.vaporPhase,),

            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            "entr_mol_ig_comp": NIST,

            "pressure_sat_comp": None,  # 修正：永久气体不计算饱和压
            "parameter_data": {
                "mw": (2.01588e-3, pyunits.kg / pyunits.mol),
                "pressure_crit": (13e5, pyunits.Pa),
                "temperature_crit": (33.18, pyunits.K),
                "cp_mol_ig_comp_coeff": {
                    'A': (33.066178, pyunits.J / (pyunits.mol * pyunits.K)),
                    'B': (-11.363417, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK)),
                    'C': (11.432816, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK ** 2)),
                    'D': (-2.772874, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK ** 3)),
                    'E': (-0.158558, (pyunits.J * pyunits.kiloK ** 2) / (pyunits.mol * pyunits.K)),
                    'F': (-9.980797, (pyunits.kiloJ / pyunits.mol)),
                    'G': (172.707974, pyunits.J / (pyunits.mol * pyunits.K)),
                    'H': (1e-10, (pyunits.kiloJ / pyunits.mol))
                },
            }
        },

        "O2": {
            "type": Component,
            "elemental_composition": {"O": 2},
            "valid_phase_types": (PhaseType.vaporPhase,),

            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            "entr_mol_ig_comp": NIST,

            "pressure_sat_comp": None,  # 永久气体不计算饱和压
            "parameter_data": {
                "mw": (31.9988e-3, pyunits.kg / pyunits.mol),
                "pressure_crit": (50.43e5, pyunits.Pa),
                "temperature_crit": (154.58, pyunits.K),
                "cp_mol_ig_comp_coeff": {
                    'A': (31.32234, pyunits.J / (pyunits.mol * pyunits.K)),
                    'B': (-20.23531, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK)),
                    'C': (57.86644, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK ** 2)),
                    'D': (-36.50624, pyunits.J / (pyunits.mol * pyunits.K * pyunits.kiloK ** 3)),
                    'E': (-0.007374, (pyunits.J * pyunits.kiloK ** 2) / (pyunits.mol * pyunits.K)),
                    'F': (-8.903471, (pyunits.kiloJ / pyunits.mol)),
                    'G': (246.7945, pyunits.J / (pyunits.mol * pyunits.K)),
                    'H': (1e-10, (pyunits.kiloJ / pyunits.mol))
                },
            }
        },

        # ============================================================
        #  KOH (liquid only)
        # ============================================================
        "KOH": {
            "type": Component,
            "valid_phase_types": (PhaseType.liquidPhase,),
            "elemental_composition": {"K": 1, "O": 1, "H": 1},

            "dens_mol_liq_comp": ConstantKOHRho,
            "enth_mol_liq_comp": ConstantKOHEnthalpy,

            "parameter_data": {
                "mw": (56.1056e-3, pyunits.kg / pyunits.mol),
            },
        },
    },

    "phases": {
        "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
        "Liq": {"type": LiquidPhase, "equation_of_state": Ideal},
    },

    # --- 相平衡全局配置 ---
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": IdealBubbleDew,

    "state_definition": FpcTP,


    "state_bounds": {
        "flow_mol_phase_comp": (1e-20, 1, 20000, pyunits.mol / pyunits.s),
        "temperature": (300, 353.15, 400, pyunits.K),
        "pressure": (1e-10, 5e5, 1e7, pyunits.Pa)
    },
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
    # "default_scaling_factors": {
    #     ("temperature", None): 1e-2,
    #     ("pressure", None): 1e-5,     # 将 Pa 缩放至 bar 量级
    #
    #     # 流量缩放 (基于预期值)
    #     ("flow_mol_phase_comp", ("Liq", "H2O")): 1e-3, # ~800 mol/s
    #     ("flow_mol_phase_comp", ("Liq", "KOH")): 1e-2,
    #     ("flow_mol_phase_comp", ("Vap", "H2")): 1e-1,  # ~13 mol/s
    #     ("flow_mol_phase_comp", ("Vap", "O2")): 1e-1,
    #     ("flow_mol_phase_comp", ("Vap", "H2O")): 1,
    #
    #     # 关键：焓缩放。目标：(流量缩放) * (焓缩放) = (功率缩放 1e-6)
    #     # 1e-3 * 1e-3 = 1e-6
    #     ("enth_mol", None): 1e-3,
    #     ("enth_mol_phase", "Liq"): 1e-3,
    #     ("enth_mol_phase", "Vap"): 1e-3,
    # }
}

