"""
Contains formuale definitions of Jacobians (in Analytical Domain)
"""
from typing import Optional, Tuple
from abc import ABC, abstractmethod
import pandas as pd
import numpy as np
from inverse_modelling_tfo.tools.s_based_intensity_datagen import get_mu_a
from tfo_sensitivity.jacobian.mu_a_equations import FullBloodJacobianMuAEqn
from tfo_sensitivity.calculate_intensity.photon_manipulation import generate_intensity_column
from .base import JacobianMuAEqn, JacobianCalculator, OperatingPoint, DxTypes


class AnalyticalJacobianCalculator(JacobianCalculator, ABC):
    """Abstract class to define how to calculate the Jacobian. Each implementation comes with its own set of
    equations to calculate different types of partial derivatives(i.e., Jacboians)"""

    def __init__(
        self,
        filtered_photon_data: pd.DataFrame,
        operating_point: OperatingPoint,
        sdd_index: Optional[int],
        base_mu_map: dict,
        dx: DxTypes,
        mu_a_eqn: Optional[JacobianMuAEqn] = None,
        normalize_derivative: bool = True,
    ) -> None:
        """
        Abstract class to define how to calculate the Jacobian(on a single detector) in a numerical way.

        Each implementation comes with its own set of equations for different types of partial derivatives
        Args:
            filtered_photon_data (pd.DataFrame): Raw photon data(should include partial paths at each layer and SDD)
            sdd_index (Optional[int]): Index of the SDD to use (Used for detector count normalization, pass None to
            skip normalization)
            base_mu_map (dict): Base mu map to use (mu_map for the non_pulsatile layers)
            dx (DxTypes): What type of partial differential to calculate
            operating_point (OperatingPoint): Operating point to use for calculating the Jacobian
            mu_a_eqn (Optional[JacobianMuAEqn], optional): Which equation to use for calculating mu_a for the fetal
            and maternal layers(Layer 1 and 4. Layer number is hardcoded). Defaults to None.
            normalize_derivative (bool, optional): Normalize the intensity derive wrt to intensity. Defaults to True.
        """
        super().__init__(
            filtered_photon_data,
            operating_point,
            sdd_index,
            base_mu_map,
            dx,
            mu_a_eqn,
        )
        # Extract the operating point values and save then individually
        self.maternal_hb = operating_point.maternal_hb
        self.maternal_sat = operating_point.maternal_sat
        self.fetal_hb = operating_point.fetal_hb
        self.fetal_sat = operating_point.fetal_sat
        self.wave_int = operating_point.wave_int
        self.mu_map = base_mu_map.copy()  # mu_map to use - will be modified by the mu_a_eqn
        if mu_a_eqn is None:
            mu_a_eqn = FullBloodJacobianMuAEqn()
        self.mu_a_eqn = mu_a_eqn
        self.normalize_derivative = normalize_derivative

        # Call initiation functions
        self._modify_mu_map()
        self.eps, self.eps_hbo, self.eps_hhb = self._calculate_eps()

    def _modify_mu_map(self) -> None:
        """
        Modify the the mu_map with the given maternal & fetal parameters
        """
        self.mu_map = self.mu_a_eqn.get_mu_map(self.base_mu_map, self.operating_point)

    def _calculate_eps(self) -> Tuple[float]:
        """
        Calculate Epsilon(extinction coefficient) for HbO and HHb
        """
        eps_hbo = get_mu_a(1.0, 1, self.wave_int)
        eps_hhb = get_mu_a(0.0, 1, self.wave_int)
        eps = get_mu_a(self.fetal_sat, 1, self.wave_int)
        return eps, eps_hbo, eps_hhb

    @abstractmethod
    def calculate_jacobian(self) -> float:
        """
        Calculate the jacobian by routing the call to the approapriate function(based on DxType)
        """

    def __str__(self) -> str:
        return f"Analytical Jacobian Calculator for {self.dx} with {self.mu_a_eqn}"


class FullBloodAnalyticalJC(AnalyticalJacobianCalculator):
    """
    Analytically calculate the Jacobian for when we assume the entire tissue is filled with blood(in terms of mu_a)
    """

    def _not_implemented(self) -> float:
        raise NotImplementedError()

    def _fetal_conc_derivative(self) -> float:
        """
        Calculate the derivative of the intensity with respect to fetal concentration
        """
        intensity_column = generate_intensity_column(self.filtered_photon_data, self.mu_map, self.sdd_index)
        pathlength_column = self.filtered_photon_data["L4 ppath"].to_numpy()
        analytical_term = -self.eps * np.dot(pathlength_column, intensity_column)
        if self.normalize_derivative:
            analytical_term /= np.sum(intensity_column)
        return analytical_term

    def _fetal_sat_derivative(self) -> float:
        intensity_column = generate_intensity_column(self.filtered_photon_data, self.mu_map, self.sdd_index)
        pathlength_column = self.filtered_photon_data["L4 ppath"].to_numpy()
        analytical_term = -(self.eps_hbo - self.eps_hhb) * self.fetal_hb * np.dot(pathlength_column, intensity_column)
        if self.normalize_derivative:
            analytical_term /= np.sum(intensity_column)
        return analytical_term

    def calculate_jacobian(self) -> float:
        # Update these values before calculating the jacobian
        self._modify_mu_map()
        self.eps, self.eps_hbo, self.eps_hhb = self._calculate_eps()
        dx_to_func_mapping = {
            "FC": self._fetal_conc_derivative,
            "FS": self._fetal_sat_derivative,
        }
        function_to_use = dx_to_func_mapping.get(self.dx, self._not_implemented)
        return function_to_use()


class PartialBloodAnalyticalJC(AnalyticalJacobianCalculator):
    """
    Analytically calculate the Jacobian for when we assume partially filled with blood(in terms of mu_a)

    Note: This assumes arterial and venous blood are the same. (This makes the equations simpler)
    Additional Args:
        arterial_volume_fraction (float): Fraction of the tissue that is filled with arterial blood. Assumes the venous
        fraction is also the same
        venous_saturation_reduction_factor (float): What fraction is the venous saturation comparted to (arterial)
        saturation
    """

    def __init__(
        self,
        filtered_photon_data: pd.DataFrame,
        operating_point: OperatingPoint,
        sdd_index: Optional[int],
        base_mu_map: dict,
        dx: DxTypes,
        mu_a_eqn: Optional[JacobianMuAEqn] = None,
        normalize_derivative: bool = True,
        arterial_volume_fraction: float = 0.1,
        # venous_volume_fraction: float = 0.1,
        venous_saturation_reduction_factor: float = 0.75,
    ) -> None:
        super().__init__(
            filtered_photon_data,
            operating_point,
            sdd_index,
            base_mu_map,
            dx,
            mu_a_eqn,
            normalize_derivative,
        )
        self.arterial_volume_fraction = arterial_volume_fraction
        # self.venous_volume_fraction = venous_volume_fraction
        self.venous_saturation_reduction_factor = venous_saturation_reduction_factor

    def _not_implemented(self) -> float:
        raise NotImplementedError()

    def _fetal_conc_derivative(self) -> float:
        """
        Calculate the derivative of the intensity with respect to fetal concentration
        """
        intensity_column = generate_intensity_column(self.filtered_photon_data, self.mu_map, self.sdd_index)
        pathlength_column = self.filtered_photon_data["L4 ppath"].to_numpy()
        eps_at_venous_sat = get_mu_a(self.fetal_sat * self.venous_saturation_reduction_factor, 1, self.wave_int)
        del_mu_del_cf = self.arterial_volume_fraction * (self.eps + eps_at_venous_sat)
        analytical_term = -del_mu_del_cf * np.dot(pathlength_column, intensity_column)
        if self.normalize_derivative:
            analytical_term /= np.sum(intensity_column)
        return analytical_term

    def _fetal_sat_derivative(self) -> float:
        intensity_column = generate_intensity_column(self.filtered_photon_data, self.mu_map, self.sdd_index)
        pathlength_column = self.filtered_photon_data["L4 ppath"].to_numpy()
        del_mu_del_sat = (
            (self.eps_hbo - self.eps_hhb)
            * self.fetal_hb
            * (self.venous_saturation_reduction_factor + 1)
            * self.arterial_volume_fraction
        )
        analytical_term = -del_mu_del_sat * np.dot(pathlength_column, intensity_column)
        if self.normalize_derivative:
            analytical_term /= np.sum(intensity_column)
        return analytical_term

    def calculate_jacobian(self) -> float:
        # Update these values before calculating the jacobian
        self._modify_mu_map()
        self.eps, self.eps_hbo, self.eps_hhb = self._calculate_eps()

        dx_to_func_mapping = {
            "FC": self._fetal_conc_derivative,
            "FS": self._fetal_sat_derivative,
        }
        function_to_use = dx_to_func_mapping.get(self.dx, self._not_implemented)
        return function_to_use()
