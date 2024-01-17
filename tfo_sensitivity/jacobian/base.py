"""
Base classes used for Jacobian calculations
"""
from dataclasses import dataclass
from abc import ABC, abstractmethod
from typing import Dict, Optional, Literal, List
import pandas as pd

# Types: M: Maternal, F: Fetal, C: Delta in concentration, S: Delta in saturation
DxTypes = Literal["MC", "MS", "FC", "FS"]


@dataclass
class OperatingPoint:
    """
    Operating Point for the Jacobian calculation
    (Based on the notion that Jacobians/partial derivatives are different at different operating points)
    Args:
        maternal_hb (float): Maternal Hb (in g/dL)
        maternal_sat (float): Maternal Saturation (0.0 to 1.0)
        fetal_hb (float): Fetal Hb (in g/dL)
        fetal_sat (float): Fetal Saturation (0.0 to 1.0)
        wave_int (int): 1 or 2 (1 for 735nm and 2 for 850nm)
    """

    maternal_hb: float
    maternal_sat: float
    fetal_hb: float
    fetal_sat: float
    wave_int: int


class JacobianMuAEqn(ABC):
    """
    Base class for calculating the mu_a for the Jacobian calculation
    """

    @abstractmethod
    def derivative_mu_map_gen(
        self, base_mu_map: Dict[int, float], operating_point: OperatingPoint, dx: DxTypes, delta: float
    ) -> List:
        """
        Creates three mu maps for calculating intensity derivatives. The first one is the base map
        modified with the properties passed in. The second on adds a postive change of delta to mu
        affected by dx and the third one adds a negative change of delta.
        """

    @abstractmethod
    def get_mu_map(self, base_mu_map: Dict[int, float], operating_point: OperatingPoint) -> Dict[int, float]:
        """
        Modifiy the mu_a of the base_mu_map based on the operating point

        Args:
            base_mu_map (Dict[int, float]): Base mu map (key: layer, value: mu_a)
            operating_point (OperatingPoint): Operating point

        Returns:
            Dict[int, float]: Modified base mu map
        """

    @abstractmethod
    def __str__(self) -> str:
        """
        A description of how the mu_a is calculated
        """


class JacobianCalculator(ABC):
    """
    Abstract class to define how to calculate the Jacobian
    """

    def __init__(
        self,
        filtered_photon_data: pd.DataFrame,
        operating_point: OperatingPoint,
        sdd_index: Optional[int],
        base_mu_map: dict,
        dx: DxTypes,
        mu_a_eqn: Optional[JacobianMuAEqn] = None,
    ) -> None:
        self.operating_point = operating_point
        self.filtered_photon_data = filtered_photon_data
        self.sdd_index = sdd_index
        self.base_mu_map = base_mu_map
        self.dx = dx
        self.mu_a_eqn = mu_a_eqn

    @abstractmethod
    def calculate_jacobian(self) -> float:
        """Calculate the jacobian using some custom formulation"""

    @abstractmethod
    def __str__(self) -> str:
        pass
