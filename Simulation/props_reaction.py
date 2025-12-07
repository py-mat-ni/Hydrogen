# rxn.py
from pyomo.environ import units as pyunits

from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock, ConcentrationForm
)
from idaes.models.properties.modular_properties.reactions.dh_rxn import \
    constant_dh_rxn
from idaes.models.properties.modular_properties.reactions.rate_forms import \
    power_law_rate
from idaes.models.properties.modular_properties.reactions.rate_constant import \
    arrhenius

# 引用你之前建立的物性包
import props

# 配置字典
config_dict = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K
    },
    "rate_reactions": {
        "water_splitting": {
            # 1. 化学计量数 (Stoichiometry)
            # 建议改为 1 mol H2O 基准，对应 z=2 电子转移，计算更直观
            "stoichiometry": {
                ("Liq", "H2O"): -1,
                ("Vap", "H2"): 1,
                ("Vap", "O2"): 0.5
            },
            # 2. 反应热 (Heat of Reaction)
            "heat_of_reaction": constant_dh_rxn,
            "parameter_data": {
                # 吸热反应，焓变为正值 (+285.83 kJ/mol)
                "dh_rxn_ref": (285.83e3, pyunits.J / pyunits.mol),

                # 3. 动力学设置 (占位符)
                # 我们会在 Unit Model 中停用这个 Arrhenius 计算
                # 并替换为法拉第定律约束
                "arrhenius_const": (1, pyunits.mol / pyunits.m ** 3 / pyunits.s),
                "energy_activation": (0, pyunits.J / pyunits.mol)
            },
            # 必须提供的形式定义 (即使我们后面会覆盖它)
            "rate_constant": arrhenius,
            "rate_form": power_law_rate,
            "concentration_form": ConcentrationForm.moleFraction,
        }
    }
}




