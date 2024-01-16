"""
Mixins to calculate tissue absorption coefficient changes used in Jacobian calculations
"""
from inverse_modelling_tfo.tools.s_based_intensity_datagen import get_mu_a
from inverse_modelling_tfo.tools.optical_properties import get_tissue_mu_a
from .base import JacobianMuAEqn


class FullBloodJacobianMuAEqn(JacobianMuAEqn):
    """
    Assumes the entire tissue is blood
    """

    def derivative_mu_map_gen(self, base_mu_map, maternal_sat, maternal_hb, fetal_sat, fetal_hb, wave_int, dx, delta):
        mu_map = base_mu_map.copy()
        mu_map[1] = get_mu_a(maternal_sat, maternal_hb, wave_int)
        mu_map[4] = get_mu_a(fetal_sat, fetal_hb, wave_int)

        mu_map1 = mu_map.copy()
        mu_map2 = mu_map.copy()
        layer_affected = 1 if "M" in dx else 4
        affected_sat = maternal_sat if "M" in dx else fetal_sat
        affected_hb = maternal_hb if "M" in dx else fetal_hb
        if "C" in dx:
            mu_map1[layer_affected] = get_mu_a(affected_sat, affected_hb + delta, wave_int)
            mu_map2[layer_affected] = get_mu_a(affected_sat, affected_hb - delta, wave_int)
        elif "S" in dx:
            mu_map1[layer_affected] = get_mu_a(affected_sat + delta, affected_hb, wave_int)
            mu_map2[layer_affected] = get_mu_a(affected_sat - delta, affected_hb, wave_int)
        else:
            raise NotImplementedError()
        return mu_map, mu_map1, mu_map2

    def __str__(self) -> str:
        return """Assumes entire tissue is blood with the same formula for mom & fetus"""


class PartialBloodJacobianMuAEqn(JacobianMuAEqn):
    """
    Assumes parts of the tissue is blood, the rest is non-blood elements. Uses the formula from
    [inverse_modelling_tfo/tools/optical_properties/get_tissue_mu_a]
    """

    def __init__(
        self,
        maternal_blood_volume_fraction: float,
        maternal_arterial_fraction: float,
        maternal_venous_saturation_reduction_factor: float,
        fetal_blood_volume_fraction: float,
        fetal_arterial_fraction: float,
        fetal_venous_saturation_reduction_factor: float,
    ) -> None:
        self.maternal_blood_volume_fraction = maternal_blood_volume_fraction
        self.fetal_blood_volume_fraction = fetal_blood_volume_fraction
        self.maternal_arterial_fraction = maternal_arterial_fraction
        self.fetal_arterial_fraction = fetal_arterial_fraction
        self.maternal_venous_saturation_reduction_factor = maternal_venous_saturation_reduction_factor
        self.fetal_venous_saturation_reduction_factor = fetal_venous_saturation_reduction_factor

    def derivative_mu_map_gen(self, base_mu_map, maternal_sat, maternal_hb, fetal_sat, fetal_hb, wave_int, dx, delta):
        mu_map = base_mu_map.copy()
        mu_map[1] = get_tissue_mu_a(
            self.maternal_blood_volume_fraction,
            maternal_hb,
            maternal_sat,
            wave_int,
            self.maternal_arterial_fraction,
            self.maternal_venous_saturation_reduction_factor,
        )
        mu_map[4] = get_tissue_mu_a(
            self.fetal_blood_volume_fraction,
            fetal_hb,
            fetal_sat,
            wave_int,
            self.fetal_arterial_fraction,
            self.fetal_venous_saturation_reduction_factor,
        )

        mu_map1 = mu_map.copy()
        mu_map2 = mu_map.copy()
        layer_affected = 1 if "M" in dx else 4
        affected_sat = maternal_sat if "M" in dx else fetal_sat
        affected_bvf = self.maternal_blood_volume_fraction if "M" in dx else self.fetal_blood_volume_fraction
        affected_avf = self.maternal_arterial_fraction if "M" in dx else self.fetal_arterial_fraction
        affected_vsrf = (
            self.maternal_venous_saturation_reduction_factor
            if "M" in dx
            else self.fetal_venous_saturation_reduction_factor
        )
        affected_hb = maternal_hb if "M" in dx else fetal_hb
        if "C" in dx:
            mu_map1[layer_affected] = get_tissue_mu_a(
                affected_bvf, affected_hb + delta, affected_sat, wave_int, affected_avf, affected_vsrf
            )
            mu_map2[layer_affected] = get_tissue_mu_a(
                affected_bvf, affected_hb - delta, affected_sat, wave_int, affected_avf, affected_vsrf
            )
        elif "S" in dx:
            mu_map1[layer_affected] = get_tissue_mu_a(
                affected_bvf, affected_hb, affected_sat + delta, wave_int, affected_avf, affected_vsrf
            )
            mu_map2[layer_affected] = get_tissue_mu_a(
                affected_bvf, affected_hb, affected_sat - delta, wave_int, affected_avf, affected_vsrf
            )
        else:
            raise NotImplementedError()
        return mu_map, mu_map1, mu_map2

    def __str__(self) -> str:
        return "Uses the formula in [inverse_modelling_tfo/tools/optical_properties/get_tissue_mu_a]"
