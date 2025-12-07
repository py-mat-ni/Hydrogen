from idaes.core import Component, VaporPhase, PhaseType
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.state_definitions import FpcTP
from idaes.models.properties.modular_properties.pure import NIST, RPP4
from pyomo.environ import units as pyunits

# --- 气相配置字典 ---
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
        #  H2O (Vapor)
        # ============================================================
        "H2O": {
            "type": Component,
            "elemental_composition": {"H": 2, "O": 1},
            "valid_phase_types": (PhaseType.vaporPhase,),  # 仅气相

            # 气相属性
            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            "entr_mol_ig_comp": NIST,
            "pressure_sat_comp": NIST,

            "parameter_data": {
                "mw": (18.0153e-3, pyunits.kg / pyunits.mol),
                "pressure_crit": (220.64e5, pyunits.Pa),
                "temperature_crit": (647.14, pyunits.K),

                # 气相热容/焓系数 (NIST)
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
                "pressure_sat_comp_coeff": {
                    'A': (4.6543, None),
                    'B': (1453.264, pyunits.K),
                    'C': (-64.848, pyunits.K),
                },
            },
        },

        # ============================================================
        #  H2 (Vapor)
        # ============================================================
        "H2": {
            "type": Component,
            "elemental_composition": {"H": 2},
            "valid_phase_types": (PhaseType.vaporPhase,),

            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            "entr_mol_ig_comp": NIST,
            "pressure_sat_comp": None,  # 永久气体

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

        # ============================================================
        #  O2 (Vapor)
        # ============================================================
        "O2": {
            "type": Component,
            "elemental_composition": {"O": 2},
            "valid_phase_types": (PhaseType.vaporPhase,),

            "enth_mol_ig_comp": NIST,
            "cp_mol_ig_comp": NIST,
            "entr_mol_ig_comp": NIST,
            "pressure_sat_comp": None,

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
    },

    "phases": {
        "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
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