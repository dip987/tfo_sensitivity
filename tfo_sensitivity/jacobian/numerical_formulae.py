"""
Contains formuale definitions of Jacobians
"""
from typing import Dict, Literal, Optional
from abc import abstractmethod, ABC
import pandas as pd
from numpy import log10
from tfo_sensitivity.calculate_intensity.photon_manipulation import generate_intensity
from tfo_sensitivity.jacobian.mu_a_equations import JacobianMuAEqn, FullBloodJacobianMuAEqn
from .base import JacobianCalculator, DxTypes, OperatingPoint


class NumericalJacobianCalculator(JacobianCalculator, ABC):
    """
    Abstract class for calculating the Jacobian(on a single detector) in a numerical way

    Args:
        filtered_photon_data (pd.DataFrame): Raw photon data(should include partial paths at each layer and SDD)
        sdd_index (Optional[int]): Index of the SDD to use (Used for detector count normalization, pass None to
        skip normalization)
        base_mu_map (dict): Base mu map to use (mu_map for the non_pulsatile layers)
        dx (DxTypes): What type of partial differential to calculate
        operating_point (OperatingPoint): Operating point to use for calculating the Jacobian
        mu_a_eqn (Optional[JacobianMuAEqn], optional): Which equation to use for calculating mu_a for the fetal
        and maternal layers(Layer 1 and 4. Layer number is hardcoded). Defaults to None.
        delta (float, optional): How much to change the values during calculating derivatives. Defaults to 0.0001.
    """

    def __init__(
        self,
        filtered_photon_data: pd.DataFrame,
        operating_point: OperatingPoint,
        sdd_index: Optional[int],
        base_mu_map: dict,
        dx: DxTypes,
        mu_a_eqn: Optional[JacobianMuAEqn] = None,
        delta: float = 0.0001,
    ) -> None:
        super().__init__(
            filtered_photon_data,
            operating_point,
            sdd_index,
            base_mu_map,
            dx,
            mu_a_eqn,
        )
        self.delta = delta
        # Extract the operating point values and save them individually
        self.maternal_hb = operating_point.maternal_hb
        self.maternal_sat = operating_point.maternal_sat
        self.fetal_hb = operating_point.fetal_hb
        self.fetal_sat = operating_point.fetal_sat
        self.wave_int = operating_point.wave_int
        # Create the mu maps for numerical derivatives
        self.mu_map = None
        self.mu_map1 = None
        self.mu_map2 = None
        if mu_a_eqn is None:
            mu_a_eqn = FullBloodJacobianMuAEqn()
        self.mu_a_eqn = mu_a_eqn
        self._derivative_mu_map_gen()

    def _derivative_mu_map_gen(self):
        """Creates three mu maps for calculating intensity derivatives. The first one is the base map
        modified with the properties passed in. The second on adds a postive change of delta to mu
        affected by dx and the third one adds a negative change of delta.
        """
        self.mu_map, self.mu_map1, self.mu_map2 = self.mu_a_eqn.derivative_mu_map_gen(
            self.base_mu_map,
            self.maternal_sat,
            self.maternal_hb,
            self.fetal_sat,
            self.fetal_hb,
            self.wave_int,
            self.dx,
            self.delta,
        )

    @abstractmethod
    def calculate_jacobian(self) -> float:
        """Calculate the jacobian using some custom formulation"""

    @abstractmethod
    def __str__(self) -> str:
        pass


class NormalizedDerivative(NumericalJacobianCalculator):
    """Calculate the partial normalized derivative from a given simulation raw data set. Uses the
    formula (I(x + delta) - I(x - delta))/I(x)/2/delta)
    """

    def calculate_jacobian(self) -> float:
        I1 = generate_intensity(self.filtered_photon_data, self.mu_map1, self.sdd_index)
        I2 = generate_intensity(self.filtered_photon_data, self.mu_map2, self.sdd_index)
        I0 = generate_intensity(self.filtered_photon_data, self.mu_map, self.sdd_index)
        return (I1 - I2) / I0 / (2 * self.delta)

    def __str__(self) -> str:
        return """Calculate the partial normalized derivative from a given simulation raw data set
        Uses the formula (I(x + delta) - I(x - delta))/I(x)/2/delta)"""


class RegularDerivative(NumericalJacobianCalculator):
    """Calculate a regular derivative with the formula
    (I(x + delta) - I(x - delta))/2/delta)
    """

    def calculate_jacobian(self) -> float:
        I1 = generate_intensity(self.filtered_photon_data, self.mu_map1, self.sdd_index)
        I2 = generate_intensity(self.filtered_photon_data, self.mu_map2, self.sdd_index)
        return (I1 - I2) / (2 * self.delta)

    def __str__(self) -> str:
        return """Calculate a regular derivative with the formula
        (I(x + delta) - I(x - delta))/2/delta)"""


class LogDerivative(NumericalJacobianCalculator):
    """Calculate a regular derivative with the formula
    (log(I(x + delta)) - log(I(x - delta)))/2/delta)
    """

    def calculate_jacobian(self) -> float:
        I1 = generate_intensity(self.filtered_photon_data, self.mu_map1, self.sdd_index)
        I2 = generate_intensity(self.filtered_photon_data, self.mu_map2, self.sdd_index)
        return (log10(I1) - log10(I2)) / (2 * self.delta)

    def __str__(self) -> str:
        return """Calculate a regular derivative with the formula
        (log(I(x + delta)) - log(I(x - delta)))/2/delta)"""


derivative_mapping: Dict[str, NumericalJacobianCalculator]
derivative_mapping = {"regular": RegularDerivative, "norm_der": NormalizedDerivative, "log": LogDerivative}
DerivativeTypes = Literal["regular", "norm_der", "log"]


def calculate_jacobian_numerical(
    derivative_format: DerivativeTypes,
    photon_data: pd.DataFrame,
    sdd_index: int,
    base_mu_map: dict,
    delta: float,
    dx: DxTypes,
    sdd_list: Optional[list],
    maternal_hb: float,
    maternal_sat: float,
    fetal_hb: float,
    fetal_sat: float,
    wave_int: int,
) -> float:
    """Convinience function to calculate the Jacobian using some given format

    Args:
        derivative_format (Literal[&#39;regular&#39;, &#39;norm_der&#39;, &#39;log&#39;]): How
        should the derivative be calculated
        photon_data (pd.DataFrame): Raw simulation file loaded from disk
        sdd_index (int): detector index to filter the data by
        base_mu_map (dict): map containing mu values for all layers. Layer 1 & 4 will be modified
        based on the given values
        delta (float): how much to change the values during calculating derivatives
        dx (Literal[&#39;MC&#39;, &#39;MS&#39;, &#39;FC&#39;, &#39;FS&#39;]): dx for partial
        derivative. M/F -> maternal or fetal layer(1 or 4), S/C -> change saturation or concetration
        sdd_list (Optional[list]): pass a list of all sdd values(sorted) to speed things up
        maternal_hb (float): current values
        maternal_sat (float): current values
        fetal_hb (float): current values
        fetal_sat (float): current values
        wave_int (int): 1 -> 735nm, 2-> 850nm

    Returns:
        float: partial derivative for dx at that SDD and for those Tissue Model Parameters
    """
    if sdd_list is None:
        sdd_list = photon_data["SDD"].unique()
    sdd = sdd_list[sdd_index]
    filtered_photon_data = photon_data[photon_data["SDD"] == sdd]
    operating_point = OperatingPoint(
        maternal_hb=maternal_hb,
        maternal_sat=maternal_sat,
        fetal_hb=fetal_hb,
        fetal_sat=fetal_sat,
        wave_int=wave_int,
    )

    jacobian_calculator = derivative_mapping[derivative_format](
        filtered_photon_data, operating_point, sdd_index, base_mu_map, delta, dx
    )
    return jacobian_calculator.calculate_jacobian()
