from pyomo.environ import units as pyunits

from idaes.models.properties.modular_properties.reactions.dh_rxn import \
    constant_dh_rxn

config_dict = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "rate_reactions": {
        # 阳极：2 OH- -> 0.5 O2 + H2O + 2e-
        "anode": {
            "stoichiometry": {
                ("Liq", "OH_ion"): -2,
                ("Vap", "O2"): 0.5,
                ("Liq", "H2O"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "parameter_data": {
                "dh_rxn_ref": (0.0, pyunits.kJ / pyunits.mol),
            },
        },
    },
}


