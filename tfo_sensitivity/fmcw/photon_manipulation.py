"""
Manipulate simulation dataframes containing raw photon data
"""
from math import ceil
from pandas import DataFrame
from numpy import linspace, digitize


def create_quantized_tof(photon_data: DataFrame, bin_count: int) -> None:
    """Create a ToF column from a give RAW photon data dataframe, inplace.
    Quantizes the ToF using the given bin count for easier plotting. You would expect and actual 
    system to have some sort of similar time resolution averaging. 

    This assumes [photon_data] has "L#_ppath" columns containing photon paths at each medium. Also, 
    the refractive index(speed of light) is constant throughout all the mediums.

    Args:
        photon_data (DataFrame): RAW photon data 
        bin_count (int): How many bins to include
    """
    def column_name_filter(column_name) -> bool:
        return ("L" in column_name) and ('ppath' in column_name)
    partial_path_columns = filter(column_name_filter, photon_data.columns)
    photon_data["Total Path"] = photon_data.loc[:,
                                                partial_path_columns].sum(axis=1)
    photon_data["ToF"] = 1/photon_data["Total Path"]
    min_tof = photon_data["ToF"].min()
    max_tof = photon_data["ToF"].max()
    quantiles = linspace(min_tof, max_tof, bin_count + 1, endpoint=True)
    photon_data["ToF"] = digitize(photon_data["ToF"], quantiles)


def create_quantized_tof_const_res(photon_data: DataFrame, time_resolution: float) -> None:
    # TODO : Debug this 
    """Create a ToF column from a give RAW photon data dataframe, inplace.
    Quantizes the ToF using the time resolution for easier plotting. So every photon within each 
    time window(of size = [time_resolution]) gets grouped together.  

    This assumes [photon_data] has "L#_ppath" columns containing photon paths at each medium. Also, 
    the refractive index(speed of light) is constant throughout all the mediums.

    Args:
        photon_data (DataFrame): RAW photon data 
        time_resolution (float): Size of each time window in mm/speed of light
    """
    def column_name_filter(column_name) -> bool:
        return ("L" in column_name) and ('ppath' in column_name)
    partial_path_columns = filter(column_name_filter, photon_data.columns)
    photon_data["Total Path"] = photon_data.loc[:,
                                                partial_path_columns].sum(axis=1)
    photon_data["ToF"] = 1/photon_data["Total Path"]
    min_tof = photon_data["ToF"].min()
    max_tof = photon_data["ToF"].max()
    required_bin_count = ceil(max_tof - min_tof/time_resolution)
    quantiles = linspace(min_tof, max_tof, required_bin_count + 1, endpoint=True)
    photon_data["ToF"] = digitize(photon_data["ToF"], quantiles)