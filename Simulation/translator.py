


from pyomo.environ import Constraint
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
)
from pyomo.common.config import ConfigBlock, ConfigValue, Bool, In


@declare_process_block_class("CompositionTranslator")
class CompositionTranslatorData(UnitModelBlockData):

    # --------------------------
    # Define CONFIG
    # --------------------------
    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare(
        "inlet_property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Inlet property package",
        ),
    )

    CONFIG.declare(
        "inlet_property_package_args",
        ConfigBlock(implicit=True)
    )

    # outlet package
    CONFIG.declare(
        "outlet_property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Outlet property package",
        ),
    )
    CONFIG.declare(
        "outlet_property_package_args",
        ConfigBlock(implicit=True)
    )

    # mapping dict
    CONFIG.declare(
        "mapping",
        ConfigValue(
            default={},
            description="(p_out, j_out) -> (p_in, j_in) mapping",
        )
    )

    # outlet state fully defined or not
    CONFIG.declare(
        "outlet_state_defined",
        ConfigValue(default=True, domain=Bool)
    )

    # ---------------------------
    # Build model
    # ---------------------------
    def build(self):
        super().build()

        # ---- Construct state blocks ----
        self.properties_in = self.config.inlet_property_package.build_state_block(
            self.flowsheet().time,
            defined_state=True,
            has_phase_equilibrium=False,
            **self.config.inlet_property_package_args
        )

        self.properties_out = self.config.outlet_property_package.build_state_block(
            self.flowsheet().time,
            defined_state=self.config.outlet_state_defined,
            has_phase_equilibrium=False,
            **self.config.outlet_property_package_args
        )

        self.add_port("inlet", self.properties_in)
        self.add_port("outlet", self.properties_out)

        mapping = self.config.mapping

        # ---------------------------
        # Link temperature & pressure
        # ---------------------------
        @self.Constraint(self.flowsheet().time)
        def eq_T(b, t):
            return b.properties_out[t].temperature == b.properties_in[t].temperature

        @self.Constraint(self.flowsheet().time)
        def eq_P(b, t):
            return b.properties_out[t].pressure == b.properties_in[t].pressure

        # ---------------------------
        # Component-by-component mapping
        # ---------------------------
        def _flow_link(b, t, p_out, j_out):
            if (p_out, j_out) in mapping:
                p_in, j_in = mapping[(p_out, j_out)]
                return (
                    b.properties_out[t].flow_mol_phase_comp[p_out, j_out]
                    ==
                    b.properties_in[t].flow_mol_phase_comp[p_in, j_in]
                )
            else:
                return b.properties_out[t].flow_mol_phase_comp[p_out, j_out] == 0

        self.eq_flow = Constraint(
            self.flowsheet().time,
            self.properties_out.phase_component_set,
            rule=_flow_link
        )
