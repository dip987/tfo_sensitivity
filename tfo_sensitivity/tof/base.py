"""
Contains the ToF class and all its corresponding functions
"""

from typing import Dict, Optional, Callable
import pandas as pd
import numpy as np
from tfo_sensitivity.calculate_intensity import generate_intensity_column

SPEED_OF_LIGHT = 3e8 / 1.4  # in m/s


class ToF:
    """
    A class to store Intensity vs. Time of Flight data and apply various operations on it.
    """

    def __init__(
        self,
        photon_data: pd.DataFrame,
        time_resolution: float,
        mu_map: Dict[int, float],
        sdd_index: Optional[int] = None,
        lower_intensity_bound: float = 1e-30,
    ):
        """
        Create a Intensity vs ToF object capable of applying different processing steps on the data.
        Args:
            photon_data (pd.DataFrame): RAW photon data (Must include an "SDD" column)
            time_resolution (float): Time resolution of the ToF (in seconds)
            mu_map (Dict[int, float]): A map containing the absorption co-eff of each layer
            sdd_index (Optional[int]): Pass the detector index to normalize with per detector count and photon count.
            Pass None to leave the intensity as is (No normalization)
            lower_intensity_bound (float): Intensity bins with values lower than this will not be stored
        """
        self.time_resolution = time_resolution
        temp_data = photon_data.copy()
        temp_data["Intensity"] = generate_intensity_column(temp_data, mu_map, sdd_index)
        temp_data["ToF_quantized"] = self._quantize_tof(self._calculate_tof(temp_data), time_resolution)
        self.photon_count_per_bin = temp_data.groupby("ToF_quantized", sort=True)["Intensity"].count()
        self.data = temp_data.groupby("ToF_quantized", sort=True)["Intensity"].sum()
        # Remove low intensity bins
        self.photon_count_per_bin = self.photon_count_per_bin.loc[self.data > lower_intensity_bound]
        self.data = self.data.loc[self.data > lower_intensity_bound]
        self.interpolation_needed = False
        self.lower_bin_index = self.data.index.values.min()
        self.upper_bin_index = self.data.index.values.max()
        self._fill_data_gaps()

    @classmethod
    def from_data(cls, data: pd.Series, time_resolution: float):
        """
        Create a ToF object from a given data series. This bypasses the initialization process and directly sets the
        internal data to a give data pandas Series
        Args:
            data (pd.Series): The data series to be used as the ToF (Intensity as the values and ToF as the index)
            time_resolution (float): Time resolution of the ToF (in seconds)
        """
        # remove index name from data - this sometimes causes naming clashes with the "ToF_quantized" column in the
        # __init__ of the ToF class (As data from another ToF would have this name for the index as well)
        data.index.name = None

        dummy_sdd = [1.0] * len(data)
        converted_pathlength = -np.log(data) * time_resolution
        new_tof = ToF(pd.DataFrame({"SDD": dummy_sdd, "L1 ppath": converted_pathlength}), time_resolution, {1: 1}, None)
        new_tof.data = data
        new_tof.lower_bin_index = new_tof.data.index.values.min()
        new_tof.upper_bin_index = new_tof.data.index.values.max()
        return new_tof

    def __truediv__(self, other):
        return self._operate(other, lambda x, y: x / y)

    def __mul__(self, other):
        return self._operate(other, lambda x, y: x * y)

    def __add__(self, other):
        return self._operate(other, lambda x, y: x + y)

    def __sub__(self, other):
        return self._operate(other, lambda x, y: x - y)

    def _operate(self, other, operation: Callable):
        if isinstance(other, ToF):
            if not self.check_operation_compatibility(other):
                raise ValueError("The two ToF objects are not compatible")
            # Keep the overlapping bins
            common_bins = self.data.index.intersection(other.data.index)
            common_data1 = self.data.loc[common_bins]
            common_data2 = other.data.loc[common_bins]
            new_data = operation(common_data1, common_data2)
            new_tof = ToF.from_data(new_data, self.time_resolution)
            return new_tof
        else:
            return operation(self.data, other)

    def __len__(self) -> int:
        return len(self.data)

    @staticmethod
    def _calculate_tof(photon_data: pd.DataFrame) -> np.ndarray:
        """Calculate the ToF column for the given RAW photon data
        Args:
            photon_data (pd.DataFrame): RAW photon data
        Returns:
            numpy array: ToF column
        """

        def column_name_filter(column_name) -> bool:
            return ("L" in column_name) and ("ppath" in column_name)

        partial_path_columns = filter(column_name_filter, photon_data.columns)
        total_path = photon_data.loc[:, partial_path_columns].sum(axis=1).to_numpy()
        # Total path is in mm -> convert to m
        tof = total_path * 1e-3 / SPEED_OF_LIGHT  # ToF in seconds
        return tof

    @staticmethod
    def _quantize_tof(tof: np.ndarray, time_resolution: float) -> np.ndarray:
        """Quantize the ToF column to the given time resolution. Returns a digit corresponding to each ToF, starting
        from 0.
        Args:
            tof (np.ndarray): ToF column
            time_resolution (float): Time resolution of the ToF (in seconds)
        Returns:
            numpy array: Quantized ToF column
        """
        return np.floor(tof / time_resolution).astype(int)

    def _fill_data_gaps(self):
        """Fill missing ToF bins with exponential, log-linear interpolated intensity"""
        lower_limit = self.data.index.values.min()
        upper_limit = self.data.index.values.max() + 1  # in the the "upper limit is not included" sense
        expected_length = upper_limit - lower_limit
        if len(self.data) < expected_length:
            self.data = self.data.reindex(range(lower_limit, upper_limit))
            self.photon_count_per_bin = self.photon_count_per_bin.reindex(range(lower_limit, upper_limit))
            self.data = np.log(self.data)
            self.photon_count_per_bin = np.log(self.photon_count_per_bin)
            self.data = self.data.interpolate(method="linear", limit_direction="both")
            self.photon_count_per_bin = self.photon_count_per_bin.interpolate(method="linear", limit_direction="both")
            self.data = np.exp(self.data)
            self.photon_count_per_bin = np.exp(self.photon_count_per_bin)
            self.interpolation_needed = True

    def check_operation_compatibility(self, other) -> bool:
        """
        Check if two ToF objects are compatible for any types of operation (Add, Div, Product, etc.)
        """
        # Check if the time resolution is the same
        if self.time_resolution != other.time_resolution:
            return False
        # Check if there is a non-zero overlap between the two ToF bins
        if self.lower_bin_index > other.upper_bin_index:
            return False
        if self.upper_bin_index < other.lower_bin_index:
            return False
        return True
