# Finding Out A Term Analogous to SNR in Simulation
SNR is an important term in maximizing sensitivity. In a noiseless setup, maximizing fetal sensitivity is analogous to always choosing the last possible, singular timebin. However, this bin would be riddled with noise in a real setup. Appending more bins to the left would bring a reduction in noise at the cost of sensitivity. The interest thus lies in optimizing both SNR & Sensitivity simultaneously. In this report, we try to investigate what term could be a proxy for SNR that can be obtained from simulation data. 


## Option 1 : Simulated Photon Count per Bin & Multiply SNR with Sensitvity
Reference to this [paper](https://opg.optica.org/oe/fulltext.cfm?uri=oe-25-23-28567&id=376623), SNR turns out to the photon count at the sample arm. Specified here as $N_S$. Simulated photon count sum over all the bins might seem like a good analogous term. However, remember in that in physics, $N \propto Intensity$. Which is not the case for our specific kind of MC simulation. Where we assume each photon has a different intensity(Not the case in real life). Which might make this a bad choice?

Regardless, here are some simulation results ![](figures/heatmap_snr1.png)
The optima occurs at [6, 127]
Adding them instead(With the SNR part normalized to 1) ![](figures/heatmap_snr2.png)
New optima at [0, 127]




## Option 2 : Intensity Sum over The Timebins & Multiply SNR with Sensitvity
Using the same logic as the section above, $N_S \propto \sum_{L_{min}}^{L_{max}} I(L)$. Using this definition of SNR, if we want to maximize $ Sensitivity \times SNR$, our equation turns into 
$$
\frac{\int I(L) R(L) dL}{\int I(L) dL} \times \int I(L)dL = \int I(L) R(L) dL
$$
The original goal of SNR was to pull the optimization window to the left. This, however, pulls it too far. From simulation data, the optimization point for this setup starts at the leftmost window and ends at some point a bit before the last timebin. Too extreme of an effect. 
