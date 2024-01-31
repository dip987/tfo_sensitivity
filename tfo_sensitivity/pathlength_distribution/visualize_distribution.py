# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tfo_sensitivity.data import load_raw
import seaborn as sns

# Plot Setup
# plt.style.use('seaborn')

# Load Sim Data
maternal_wall_thickness, uterus_thickness, wave_int = 12, 5, 2
raw_sim_data_path = load_raw(maternal_wall_thickness, uterus_thickness, wave_int)
sim_data = pd.read_pickle(raw_sim_data_path)
# Create SDD column!
sim_data['SDD'] = sim_data['X'] - 100
all_sdd = sim_data['SDD'].unique()

all_sdd = sim_data['SDD'].unique()
all_sdd.sort()
current_sdd = all_sdd[1]

sim_data = sim_data[sim_data['SDD'] == current_sdd]
sim_data["sum_L"] = sim_data["L1 ppath"] + sim_data["L2 ppath"] + sim_data["L3 ppath"]
sim_data.describe()

# Create a 2D histogram of the sum_L4 and L4 ppath columns
bin_edges_x = np.linspace(0, 2000, 200)
bin_edges_y = np.linspace(0, 3000, 200)
hist_2d, x_edges, y_edges = np.histogram2d(sim_data["L1 ppath"], sim_data["L2 ppath"], bins=[bin_edges_x, bin_edges_y],
                                           density=True)
# hist_2d, x_edges, y_edges = np.histogram2d(sim_data["sum_L"], sim_data["L4 ppath"], bins=400)

# Plot the 2D histogram
fig, ax = plt.subplots()
EPSILON = 1e-12
hist_2d[hist_2d == 0.0] = EPSILON
plt.imshow(np.log10(hist_2d.T), interpolation='nearest', origin='lower',
           extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]], aspect='auto', cmap='summer')
plt.colorbar()
plt.title("Partial Path Distribution Density(Log)")
plt.xlabel('Sum L1, L2, L3 ppath')
plt.ylabel("L4 ppath")
plt.show()
