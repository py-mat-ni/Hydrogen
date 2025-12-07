from pyomo.environ import (
    ConcreteModel,
    Var,
    Constraint,
    TransformationFactory,
    value,
    Block,
    Expression
)

from idaes.core import FlowsheetBlock
from pyomo.network import Arc
from copy import deepcopy

import idaes.logger as idaeslog
from idaes.core.solvers import get_solver

import props as props
import props_reaction as props_reaction
import rxn_cat as rxn_cat
import rxn_an as rxn_an

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
)
from idaes.models.properties.modular_properties.state_definitions import FpcTP
from idaes.models.unit_models import (
    Feed,
    Mixer,
    Separator,
    SplittingType,
    StoichiometricReactor,
    Flash,
    Product,
    Heater,
    PressureChanger,
)
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom
import time
from idaes.core.util.model_serializer import to_json, from_json



# Build model
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

# Property and reaction packages
m.fs.props = GenericParameterBlock(**props.configuration)

# 创建反应包
m.fs.rxn_cat = GenericReactionParameterBlock(
    property_package=m.fs.props,
    **rxn_cat.config_dict,
)
m.fs.rxn_an = GenericReactionParameterBlock(
    property_package=m.fs.props,
    **rxn_an.config_dict)



m.fs.Mix_C = Mixer(
    property_package=m.fs.props,
    inlet_list=["makeup", "recycle"],
    momentum_mixing_type=MomentumMixingType.equality
)
m.fs.Cathode = StoichiometricReactor(
    property_package=m.fs.props,
    reaction_package=m.fs.rxn_cat,
    has_heat_of_reaction=True,
    has_heat_transfer=True,
    has_pressure_change=False,
)
m.fs.Flash_C = Flash(
    property_package=m.fs.props,
    has_heat_transfer=False,
    has_pressure_change=True
)
m.fs.mix1 = Mixer(
    property_package=m.fs.props,
    inlet_list=["one", "two"],
    momentum_mixing_type=MomentumMixingType.equality
    )
m.fs.Sep_C = Separator(
    property_package=m.fs.props,
    outlet_list=["one", "two"],
    split_basis=SplittingType.phaseComponentFlow,
    )
m.fs.Sep_C1 = Separator(
    property_package=m.fs.props,
    outlet_list=["to_h2", "to_mix"],
    split_basis=SplittingType.phaseComponentFlow,
    )
m.fs.Sep_C_mem = Separator(
    property_package=m.fs.props,
    outlet_list=["to_mix", "to_anode"],
    split_basis=SplittingType.phaseComponentFlow,
)
m.fs.HE_H2 = Heater(
    property_package=m.fs.props, has_pressure_change=False
)
m.fs.Mix_H2_wash = Mixer(
    property_package=m.fs.props,
    inlet_list=["h2_gas", "koh_feed"],
    momentum_mixing_type=MomentumMixingType.equality,
)
m.fs.Flash_H2_wash = Flash(
    property_package=m.fs.props, has_heat_transfer=False, has_pressure_change=True
)

m.fs.Mix_KOH_C = Mixer(
    property_package=m.fs.props,
    inlet_list=["from_Sep", "from_wash", "from_OH"],
    momentum_mixing_type=MomentumMixingType.equality,
)
m.fs.Split_ab = Separator(
    property_package=m.fs.props,
    outlet_list=["a", "b"],
)
# Anode loop: Mixer (incl. from membrane) -> StoichReactor -> Flash -> recycle
m.fs.Mix_A = Mixer(
    property_package=m.fs.props,
    inlet_list=["from_cathode", "recycle"],
    momentum_mixing_type=MomentumMixingType.equality,
)


m.fs.Anode = StoichiometricReactor(
    property_package=m.fs.props,
    reaction_package=m.fs.rxn_an,
    has_heat_of_reaction=True,
    has_heat_transfer=True,
    has_pressure_change=False,
)
m.fs.Flash_A = Flash(
    property_package=m.fs.props, has_heat_transfer=False, has_pressure_change=True
)
m.fs.Sep_A = Separator(
    property_package=m.fs.props,
    outlet_list=["one","two"],
    split_basis=SplittingType.phaseComponentFlow,
    )
m.fs.mix2 = Mixer(
    property_package=m.fs.props,
    inlet_list=["one", "two"],
    momentum_mixing_type=MomentumMixingType.equality
    )
m.fs.Sep_A1 = Separator(
    property_package=m.fs.props,
    outlet_list=["to_O2","to_mix"],
    split_basis=SplittingType.phaseComponentFlow,
    )
m.fs.HE_O2 = Heater(
    property_package=m.fs.props, has_pressure_change=False
)
m.fs.Mix_O2_wash = Mixer(
    property_package=m.fs.props,
    inlet_list=["o2_gas", "koh_feed"],
    momentum_mixing_type=MomentumMixingType.equality,
)
m.fs.Flash_O2_wash = Flash(
    property_package=m.fs.props, has_heat_transfer=False, has_pressure_change=True
)
m.fs.Mix_KOH_A = Mixer(
    property_package=m.fs.props,
    inlet_list=["from_sep", "from_wash"],
    momentum_mixing_type=MomentumMixingType.equality,
)
m.fs.Split_cd = Separator(
    property_package=m.fs.props,
    outlet_list=["c", "d"],
)
m.fs.Mix_ac = Mixer(
    property_package=m.fs.props,
    inlet_list=["a", "c"],
    momentum_mixing_type=MomentumMixingType.equality,
)
m.fs.HE_KOH1 = Heater(property_package=m.fs.props, has_pressure_change=False)
m.fs.Press_KOH_C = PressureChanger(
    property_package=m.fs.props,
    thermodynamic_assumption=ThermodynamicAssumption.isothermal,
    compressor=True,
)
m.fs.Mix_bd = Mixer(
    property_package=m.fs.props,
    inlet_list=["b", "d"],
    momentum_mixing_type=MomentumMixingType.equality,
)
m.fs.HE_KOH2 = Heater(property_package=m.fs.props, has_pressure_change=False)
m.fs.Press_KOH_A = PressureChanger(
    property_package=m.fs.props,
    thermodynamic_assumption=ThermodynamicAssumption.isothermal,
    compressor=True,
)
# Feeds and products
m.fs.Catholyte_Feed = Feed(property_package=m.fs.props)
m.fs.Anolyte_Feed = Feed(property_package=m.fs.props)
m.fs.H2_product = Product(property_package=m.fs.props)
m.fs.O2_product = Product(property_package=m.fs.props)

# Connectivity
# Cathode
m.fs.a01 = Arc(source=m.fs.Mix_C.outlet, destination=m.fs.Cathode.inlet)
m.fs.a02 = Arc(source=m.fs.Cathode.outlet, destination=m.fs.Flash_C.inlet)

m.fs.a001 = Arc(source=m.fs.Flash_C.vap_outlet, destination=m.fs.mix1.one)
m.fs.a002 = Arc(source=m.fs.Flash_C.liq_outlet, destination=m.fs.mix1.two)

m.fs.a003 = Arc(source=m.fs.mix1.outlet, destination=m.fs.Sep_C.inlet)


m.fs.a03 = Arc(source=m.fs.Sep_C.one, destination=m.fs.Sep_C1.inlet)
m.fs.a04 = Arc(source=m.fs.Sep_C.two, destination=m.fs.Sep_C_mem.inlet)

m.fs.a031 = Arc(source=m.fs.Sep_C1.to_h2, destination=m.fs.HE_H2.inlet)
m.fs.a041 = Arc(source=m.fs.Sep_C1.to_mix, destination=m.fs.Mix_KOH_C.from_Sep)

m.fs.a03h = Arc(source=m.fs.HE_H2.outlet, destination=m.fs.Mix_H2_wash.h2_gas)
m.fs.a03f = Arc(source=m.fs.Mix_H2_wash.outlet, destination=m.fs.Flash_H2_wash.inlet)
m.fs.a03p = Arc(source=m.fs.Flash_H2_wash.vap_outlet, destination=m.fs.H2_product.inlet)




m.fs.a05m = Arc(source=m.fs.Sep_C_mem.to_anode, destination=m.fs.Mix_A.from_cathode)
m.fs.a051 = Arc(source=m.fs.Sep_C_mem.to_mix, destination=m.fs.Mix_KOH_C.from_OH)


m.fs.a06w = Arc(source=m.fs.Catholyte_Feed.outlet, destination=m.fs.Mix_H2_wash.koh_feed)
m.fs.a06l = Arc(source=m.fs.Flash_H2_wash.liq_outlet, destination=m.fs.Mix_KOH_C.from_wash)
m.fs.a06s = Arc(source=m.fs.Mix_KOH_C.outlet, destination=m.fs.Split_ab.inlet)
m.fs.a06a = Arc(source=m.fs.Split_ab.a, destination=m.fs.Mix_ac.a)
m.fs.a06b = Arc(source=m.fs.Split_ab.b, destination=m.fs.Mix_bd.b)

# Anode
m.fs.a11 = Arc(source=m.fs.Mix_A.outlet, destination=m.fs.Anode.inlet)
m.fs.a12 = Arc(source=m.fs.Anode.outlet, destination=m.fs.Flash_A.inlet)

m.fs.a12a = Arc(source=m.fs.Flash_A.vap_outlet, destination=m.fs.Sep_A.inlet)
m.fs.a12b = Arc(source=m.fs.Flash_A.liq_outlet, destination=m.fs.Sep_A.inlet)

m.fs.a112 = Arc(source=m.fs.Sep_A.one, destination=m.fs.mix2.one)
m.fs.a122 = Arc(source=m.fs.Sep_A.two, destination=m.fs.mix2.two)

m.fs.a132 = Arc(source=m.fs.mix2.outlet, destination=m.fs.Sep_A1.inlet)



m.fs.a13 = Arc(source=m.fs.Sep_A1.to_O2, destination=m.fs.HE_O2.inlet)
m.fs.a13h = Arc(source=m.fs.HE_O2.outlet, destination=m.fs.Mix_O2_wash.o2_gas)
m.fs.a13f = Arc(source=m.fs.Mix_O2_wash.outlet, destination=m.fs.Flash_O2_wash.inlet)
m.fs.a13p = Arc(source=m.fs.Flash_O2_wash.vap_outlet, destination=m.fs.O2_product.inlet)

m.fs.a14 = Arc(source=m.fs.Sep_A1.to_mix, destination=m.fs.Mix_KOH_A.from_sep)
m.fs.a15w = Arc(source=m.fs.Anolyte_Feed.outlet, destination=m.fs.Mix_O2_wash.koh_feed)
m.fs.a15l = Arc(source=m.fs.Flash_O2_wash.liq_outlet, destination=m.fs.Mix_KOH_A.from_wash)
m.fs.a15s = Arc(source=m.fs.Mix_KOH_A.outlet, destination=m.fs.Split_cd.inlet)
m.fs.a15c = Arc(source=m.fs.Split_cd.c, destination=m.fs.Mix_ac.c)
m.fs.a15d = Arc(source=m.fs.Split_cd.d, destination=m.fs.Mix_bd.d)

m.fs.r01 = Arc(source=m.fs.Mix_ac.outlet, destination=m.fs.HE_KOH1.inlet)
m.fs.r02 = Arc(source=m.fs.HE_KOH1.outlet, destination=m.fs.Press_KOH_C.inlet)
m.fs.r03 = Arc(source=m.fs.Press_KOH_C.outlet, destination=m.fs.Mix_C.recycle)
m.fs.r11 = Arc(source=m.fs.Mix_bd.outlet, destination=m.fs.HE_KOH2.inlet)
m.fs.r12 = Arc(source=m.fs.HE_KOH2.outlet, destination=m.fs.Press_KOH_A.inlet)
m.fs.r13 = Arc(source=m.fs.Press_KOH_A.outlet, destination=m.fs.Mix_A.recycle)


# Expand arcs
TransformationFactory("network.expand_arcs").apply_to(m)
print("Initial DOF =", degrees_of_freedom(m))

m.fs.Catholyte_Feed.outlet.temperature.fix(298)
m.fs.Catholyte_Feed.outlet.pressure.fix(1.5e5)
m.fs.Catholyte_Feed.outlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(8 - 3e-7)
m.fs.Catholyte_Feed.outlet.flow_mol_phase_comp[0, "Liq", "K_ion"].fix(1)
m.fs.Catholyte_Feed.outlet.flow_mol_phase_comp[0, "Liq", "OH_ion"].fix(1)
m.fs.Catholyte_Feed.outlet.flow_mol_phase_comp[0, "Vap", "H2O"].fix(1e-7)
m.fs.Catholyte_Feed.outlet.flow_mol_phase_comp[0, "Vap", "H2"].fix(1e-7)
m.fs.Catholyte_Feed.outlet.flow_mol_phase_comp[0, "Vap", "O2"].fix(1e-7)

m.fs.Anolyte_Feed.outlet.temperature.fix(298)
m.fs.Anolyte_Feed.outlet.pressure.fix(1.5e5)
m.fs.Anolyte_Feed.outlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(8 - 3e-7)
m.fs.Anolyte_Feed.outlet.flow_mol_phase_comp[0, "Liq", "K_ion"].fix(1)
m.fs.Anolyte_Feed.outlet.flow_mol_phase_comp[0, "Liq", "OH_ion"].fix(1)
m.fs.Anolyte_Feed.outlet.flow_mol_phase_comp[0, "Vap", "H2O"].fix(1e-7)
m.fs.Anolyte_Feed.outlet.flow_mol_phase_comp[0, "Vap", "H2"].fix(1e-7)
m.fs.Anolyte_Feed.outlet.flow_mol_phase_comp[0, "Vap", "O2"].fix(1e-7)

m.fs.Cathode.heat_duty.fix(0.0)
m.fs.Anode.heat_duty.fix(0.0)
m.fs.HE_H2.heat_duty.fix(-3e2)
m.fs.HE_O2.heat_duty.fix(-3e2)
m.fs.HE_KOH1.heat_duty.fix(-3e2)
m.fs.HE_KOH2.heat_duty.fix(-3e2)
m.fs.Press_KOH_C.deltaP.fix(2e5)
m.fs.Press_KOH_A.deltaP.fix(2e5)
m.fs.Split_ab.split_fraction[0, "a"].fix(0.5)
m.fs.Split_cd.split_fraction[0, "c"].fix(0.5)

def make_isothermal_isobaric_flash(flash, tag):
    setattr(m.fs, f"{tag}_T_liq_eq",
        Constraint(expr=flash.liq_outlet.temperature[0] == flash.inlet.temperature[0]))
# 对你模型里的各个 Flash 应用
make_isothermal_isobaric_flash(m.fs.Flash_C,        "Flash_C")
make_isothermal_isobaric_flash(m.fs.Flash_A,        "Flash_A")
make_isothermal_isobaric_flash(m.fs.Flash_H2_wash,  "Flash_H2_wash")
make_isothermal_isobaric_flash(m.fs.Flash_O2_wash,  "Flash_O2_wash")


m.fs.I = Var(initialize=10.0)       # A
m.fs.F_const = 96485.0              # C/mol
m.fs.cathode_faraday = Constraint(
    expr=m.fs.Cathode.control_volume.rate_reaction_extent[0, "cathode"]
         == m.fs.I / (2 * m.fs.F_const)
)
m.fs.anode_faraday = Constraint(
    expr=m.fs.Anode.control_volume.rate_reaction_extent[0, "anode"]
         == m.fs.I / (4 * m.fs.F_const)
)

# 氢气、氧气只存在气相
for j in ["H2", "O2", "H2O"]:
    m.fs.Sep_C_mem.split_fraction[0, "to_anode", "Vap", j].fix(1e-7)
for j in ["K_ion", "H2O"]:
    m.fs.Sep_C_mem.split_fraction[0, "to_anode", "Liq", j].fix(1e-7)

# 阴极侧：定义穿膜的 OH⁻ 摩尔通量
m.fs.N_OH = Expression(expr=m.fs.I / m.fs.F_const)
def _oh_split_rule(b):
    eps = 1e-12  # 防止除零
    return m.fs.Sep_C_mem.split_fraction[0, "to_anode", "Liq", "OH_ion"] == \
        m.fs.N_OH / (m.fs.Sep_C_mem.inlet.flow_mol_phase_comp[0, "Liq", "OH_ion"] + eps)
m.fs.oh_split_def = Constraint(rule=_oh_split_rule)



print("DOF before solve =", degrees_of_freedom(m))
# to_json(m, fname="AWE_flowsheet_state.json")
# m.fs.visualize("AWE_flowsheet")
# print("Cathode DOF =", degrees_of_freedom(m.fs.Cathode))
# print("Flash_C DOF =", degrees_of_freedom(m.fs.Flash_C))
# print("Sep_C DOF =", degrees_of_freedom(m.fs.Sep_C))
# print("Sep_C1 DOF =", degrees_of_freedom(m.fs.Sep_C1))
# print("Sep_C_mem DOF =", degrees_of_freedom(m.fs.Sep_C_mem))
# print("Anode DOF =", degrees_of_freedom(m.fs.Anode))
# print("Flash_A DOF =", degrees_of_freedom(m.fs.Flash_A))
# print("Sep_A DOF =", degrees_of_freedom(m.fs.Sep_A))
# print("Sep_A1 DOF =", degrees_of_freedom(m.fs.Sep_A1))
# print("HE_H2 DOF =", degrees_of_freedom(m.fs.HE_H2))
# print("HE_O2 DOF =", degrees_of_freedom(m.fs.HE_O2))
# print("HE_KOH1 DOF =", degrees_of_freedom(m.fs.HE_KOH1))
# print("HE_KOH2 DOF =", degrees_of_freedom(m.fs.HE_KOH2))
# print("Press_KOH_C DOF =", degrees_of_freedom(m.fs.Press_KOH_C))
# print("Press_KOH_A DOF =", degrees_of_freedom(m.fs.Press_KOH_A))
# print("Split_ab DOF =", degrees_of_freedom(m.fs.Split_ab))
# print("Split_cd DOF =", degrees_of_freedom(m.fs.Split_cd))
# print("Mix_ac DOF =", degrees_of_freedom(m.fs.Mix_ac))
# print("Mix_bd DOF =", degrees_of_freedom(m.fs.Mix_bd))
# print("Mix_C DOF =", degrees_of_freedom(m.fs.Mix_C))
# print("Mix_A DOF =", degrees_of_freedom(m.fs.Mix_A))
# print("Mix_H2_wash DOF =", degrees_of_freedom(m.fs.Mix_H2_wash))
# print("Mix_O2_wash DOF =", degrees_of_freedom(m.fs.Mix_O2_wash))
# print("Mix_KOH_C DOF =", degrees_of_freedom(m.fs.Mix_KOH_C))
# print("Mix_KOH_A DOF =", degrees_of_freedom(m.fs.Mix_KOH_A))
# print("Flash_H2_wash DOF =", degrees_of_freedom(m.fs.Flash_H2_wash))
# print("Flash_O2_wash DOF =", degrees_of_freedom(m.fs.Flash_O2_wash))
# print("=== Unfixed variables in Sep_C ===")
# for v in m.fs.Sep_C.component_data_objects(Var, descend_into=True):
#     if not v.fixed:
#         print(v.name)


# while True:
#     time.sleep(1)

# New units initialization (order respecting connectivity)
m.fs.HE_H2.initialize(outlvl=idaeslog.INFO_LOW)
m.fs.Mix_H2_wash.initialize(outlvl=idaeslog.INFO_LOW)
m.fs.Flash_H2_wash.initialize(outlvl=idaeslog.INFO_LOW)
m.fs.Mix_KOH_C.initialize(outlvl=idaeslog.INFO_LOW)
m.fs.Split_ab.initialize(outlvl=idaeslog.INFO_LOW)
m.fs.HE_O2.initialize(outlvl=idaeslog.INFO_LOW)
m.fs.Mix_O2_wash.initialize(outlvl=idaeslog.INFO_LOW)
m.fs.Flash_O2_wash.initialize(outlvl=idaeslog.INFO_LOW)
m.fs.Mix_KOH_A.initialize(outlvl=idaeslog.INFO_LOW)
m.fs.Split_cd.initialize(outlvl=idaeslog.INFO_LOW)
m.fs.Mix_ac.initialize(outlvl=idaeslog.INFO_LOW)
m.fs.HE_KOH1.initialize(outlvl=idaeslog.INFO_LOW)
m.fs.Press_KOH_C.initialize(outlvl=idaeslog.INFO_LOW)
m.fs.Mix_bd.initialize(outlvl=idaeslog.INFO_LOW)
m.fs.HE_KOH2.initialize(outlvl=idaeslog.INFO_LOW)
m.fs.Press_KOH_A.initialize(outlvl=idaeslog.INFO_LOW)

solver = get_solver()
solver.options["max_iter"] = 5000

# Optionally ramp current to target (example 100 A)
# m.fs.I.fix(100.0)

results = solver.solve(m, tee=True)
print("Solve status:", results.solver.status)

print("H2 flow (mol/s):", value(m.fs.H2_product.flow_mol_phase_comp[0, "Vap", "H2"]))
print("O2 flow (mol/s):", value(m.fs.O2_product.flow_mol_phase_comp[0, "Vap", "O2"]))
print("Membrane OH flux (mol/s):", value(m.fs.N_OH))
to_json(m, fname="AWE_flowsheet_state.json")
m.fs.visualize("AWE_flowsheet")
while True:
    time.sleep(1)
