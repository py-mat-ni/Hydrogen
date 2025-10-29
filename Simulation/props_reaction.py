
from pyomo.environ import units as pyunits

from idaes.models.properties.modular_properties.base.generic_reaction import (
        ConcentrationForm)
from idaes.models.properties.modular_properties.reactions.dh_rxn import \
    constant_dh_rxn
from idaes.models.properties.modular_properties.reactions.rate_forms import \
    power_law_rate
from idaes.models.properties.modular_properties.reactions.rate_constant import \
    arrhenius

config_dict = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K
    },
    "rate_reactions": {
        "H2comb1": {
            "stoichiometry": {
                ("Liq", "H2O"): -2,
                ("Vap", "H2"): 2,
                ("Vap", "O2"): 1
            },
            "heat_of_reaction": constant_dh_rxn,
            "rate_constant": arrhenius,  # gibbs_energy,
            "rate_form": power_law_rate,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                # ToDo Check values
                "dh_rxn_ref": (285.83, pyunits.kiloJ / pyunits.mol),  # Scheepers2020 286
                "arrhenius_const": (0, pyunits.mol / pyunits.m ** 3 / pyunits.s),  # for arrhenius
                "energy_activation": (0, pyunits.J / pyunits.mol),  # Lettenmeier2016 40-52 kJ/mol
                # "ds_rxn_ref": (0.1644, pyunits.kiloJ/pyunits.mol/pyunits.K),  # S=(H-G)/T=(286-237)/298 Scheepers2020
                # "temperature_eq_ref": (298.15, pyunits.K)  # (273.15 + 80, pyunits.K)
            }
        }
    },
}



