"Helper functions to calculate intensity from RAW simulation files"
from typing import Optional
from pandas import DataFrame
import numpy as np
from inverse_modelling_tfo.data import EQUIDISTANCE_DETECTOR_COUNT, EQUIDISTANCE_DETECTOR_PHOTON_COUNT


def generate_intensity(photon_data: DataFrame, mu_map: dict, sdd_index: Optional[int]) -> float:
    """Calculates the combined intensity from all of photons present in the Dataframe
    (This assumes all the proper mus are in the mu_map)
    Args:
        photon_data (DataFrame): RAW photon data
        mu_map (dict): a map containing the absorption co-eff of each layer
        sdd_index (Optional[int]): Pass the detector index to normalize the per detector count and
        photon count. Pass None to leave the intensity as is (No normalization)
    """
    intensity = generate_intensity_column(photon_data, mu_map, sdd_index)
    intensity = np.sum(intensity)
    return intensity


def generate_intensity_column(photon_data: DataFrame, mu_map: dict,
                              sdd_index: Optional[int]=None) -> np.ndarray:
    """Calculates intensity of each photon present in the Dataframe and returns it as a long vector
    (This assumes all the proper mu's are in the mu_map)
    Args:
        photon_data (DataFrame): RAW photon data
        mu_map (dict): a map containing the absorption co-eff of each layer
        sdd_index (Optional[int]): Pass the detector index to normalize the per detector count and
        photon count. Pass None to leave the intensity as is (No normalization)
    
    Returns:
        (np.ndarray): Intensity per photon
    """
    def name_filer(column_name) -> bool:
        return ("L" in column_name) and ('ppath' in column_name)
    partial_path_columns = list(filter(name_filer, photon_data.columns))
    assert len(partial_path_columns) == len(
        mu_map), "mu_map & partial path data length mismatch"

    # Sort sets
    mu_map = dict(sorted(mu_map.items()))
    partial_path_columns.sort()

    intensity = photon_data[partial_path_columns].to_numpy()
    for layer_number, mu in mu_map.items():
        intensity[:, layer_number -
                  1] = np.exp(-mu * intensity[:, layer_number - 1])
    intensity = np.prod(intensity, axis=1)
    if sdd_index is not None:
        intensity /= EQUIDISTANCE_DETECTOR_COUNT[sdd_index]
        intensity /= EQUIDISTANCE_DETECTOR_PHOTON_COUNT
    return intensity
