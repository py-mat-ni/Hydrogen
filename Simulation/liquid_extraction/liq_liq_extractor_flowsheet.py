#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
    The below is an implementation of a flowsheet for liquid liquid extractor. 
    The unit model uses two property packages for two liquid fluids
"""
import pyomo.environ as pyo
from idaes.core.solvers import get_solver
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)
from idaes.core import FlowsheetBlock
# from idaes_examples.notebooks.docs.unit_models.custom_unit_models.liquid_extraction.organic_property import (
#     OrgPhase,
# )
# from idaes_examples.notebooks.docs.unit_models.custom_unit_models.liquid_extraction.aqueous_property import (
#     AqPhase,
# )
# from idaes_examples.notebooks.docs.unit_models.custom_unit_models.liquid_extraction.liquid_liquid_extractor import (
#     LiqExtraction,
# )
from Simulation.liquid_extraction.organic_property import OrgPhase
from Simulation.liquid_extraction.aqueous_property import AqPhase
from Simulation.liquid_extraction.liquid_liquid_extractor import LiqExtraction
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.environ import Var
def build_model(m):
    m.fs.org_properties = OrgPhase()
    m.fs.aq_properties = AqPhase()

    m.fs.lex = LiqExtraction(
        dynamic=False,
        has_pressure_change=False,
        organic_property_package=m.fs.org_properties,
        aqueous_property_package=m.fs.aq_properties,
    )


def fix_initial_state(m):
    print("Initial DOF =", degrees_of_freedom(m))
    m.fs.lex.organic_inlet.flow_vol.fix(100 * pyo.units.L / pyo.units.hour)
    print("Initial DOF =", degrees_of_freedom(m))
    m.fs.lex.organic_inlet.temperature.fix(300 * pyo.units.K)
    print("Initial DOF =", degrees_of_freedom(m))
    m.fs.lex.organic_inlet.pressure.fix(1 * pyo.units.atm)
    print("Initial DOF =", degrees_of_freedom(m))
    m.fs.lex.organic_inlet.conc_mass_comp[0, "NaCl"].fix(
        1e-6 * pyo.units.g / pyo.units.L
    )
    m.fs.lex.organic_inlet.conc_mass_comp[0, "KNO3"].fix(
        1e-6 * pyo.units.g / pyo.units.L
    )
    m.fs.lex.organic_inlet.conc_mass_comp[0, "CaSO4"].fix(
        1e-6 * pyo.units.g / pyo.units.L
    )
    print("Initial DOF =", degrees_of_freedom(m))
    m.fs.lex.aqueous_inlet.flow_vol.fix(100 * pyo.units.L / pyo.units.hour)

    m.fs.lex.aqueous_inlet.temperature.fix(300 * pyo.units.K)

    m.fs.lex.aqueous_inlet.pressure.fix(1 * pyo.units.atm)
    print("Initial DOF =", degrees_of_freedom(m))
    m.fs.lex.aqueous_inlet.conc_mass_comp[0, "NaCl"].fix(
        1 * pyo.units.g / pyo.units.L
    )
    m.fs.lex.aqueous_inlet.conc_mass_comp[0, "KNO3"].fix(
        1 * pyo.units.g / pyo.units.L
    )
    m.fs.lex.aqueous_inlet.conc_mass_comp[0, "CaSO4"].fix(
        1 * pyo.units.g / pyo.units.L
    )
    print("Initial DOF =", degrees_of_freedom(m))
    # The below changes are done because the same property package is used to demonstrate diagnostics
    # The diagnostics section is illustrated using the below example.
    m.fs.org_properties.diffusion_factor["NaCl"] = (
        m.fs.org_properties.diffusion_factor["NaCl"] / 100
    )
    m.fs.org_properties.diffusion_factor["KNO3"] = (
        m.fs.org_properties.diffusion_factor["KNO3"] / 100
    )
    m.fs.org_properties.diffusion_factor["CaSO4"] = (
        m.fs.org_properties.diffusion_factor["CaSO4"] / 100
    )

    m.fs.lex.organic_phase.properties_in[0.0].pressure.setlb(0.5)
    m.fs.lex.organic_phase.properties_out[0.0].pressure.setlb(0.5)


def initialize_model(m):
    initializer = BlockTriangularizationInitializer()
    initializer.initialize(m.fs.lex)


def solve_model(m):
    solver = get_solver()
    results = solver.solve(m, tee=True)


def main():
    m = pyo.ConcreteModel(name="NGFC no CCS")
    m.fs = FlowsheetBlock(dynamic=False)
    build_model(m)
    fix_initial_state(m)
    initialize_model(m)
    solve_model(m)
    import os
    import time
    if os.path.exists("Hy.json"):
        os.remove("Hy.json")
    m.fs.visualize("Hy")
    for v in m.component_data_objects(Var, descend_into=True):
        if (not v.fixed) and (v.is_potentially_variable()):
            # 可以再加点过滤，只看 washer 里的：
            if "fs.washer" in v.name:
                print(v.name)
    while True:
        time.sleep(1)


if __name__ == "__main__":
    main()
