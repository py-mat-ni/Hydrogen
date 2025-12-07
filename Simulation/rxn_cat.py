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
        # 阴极：2 H2O + 2e- -> H2 + 2 OH-
        "cathode": {
            "stoichiometry": {
                ("Liq", "H2O"): -2,
                ("Vap", "H2"): 1,
                ("Liq", "OH_ion"): 2,
            },
            "heat_of_reaction": constant_dh_rxn,
            "parameter_data": {
                # 将反应热置 0，由电功/换热器承担能量收支
                "dh_rxn_ref": (0.0, pyunits.kJ / pyunits.mol),
            },
        },
    },
}



