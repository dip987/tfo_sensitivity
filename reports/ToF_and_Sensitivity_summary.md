# Relating Sensitivity to ToF Filters
## Defining Sensitivity 
We define sensitivity as the relative partial derivative of the intensity signal wrt any of the fetal TMP we care about(So either saturation or [Hb])
```math
Sensitivity = \frac {\delta I} {\delta TMP_{fetal}} / I 
```

## ToF Measurements & ToF Filtering
A regular ToF measurement would give us an intensity vs. time of flight for the photons, $I(T)$. Assuming a constant speed of light through every tissue layer, we can replace T with the total photon pathlength, L. This allows us to reparameterize $I(T)$ and break it down using Beer-Lambert's law
```math
I_{ToF}(T) = I_{ToF}(L) = \int_{L_1} \int_{L_2} .. \int_{L_M} \mathfrak{1}_{(\sum_jL_j = L)} exp(-\sum_j \mu_{a,j}L_j)p(L_1, L_2, .., L_M)
```
Where, j denotes each tissue layer in the model, and $\mu_{a,j}$ is its corresponding absorption coefficeint. $L_j$ denotes each photon's pathlength through only the j-th layer. $p(L_1, L_2, .., L_M)$ is the joint probability distribution over the partial optical pathlengths through each tissue layer. This distribution would mainly rely on the tissue geometry, how the layers are organized with the tissue, and the scattering co-efficient($\mu_s$) for each layer. The term $\mathfrak{1}_{(\sum_jL_j = L)}$ is 1 when the condition $\sum_jL_j = L$ is met, otherwise 0. 

For a simplified visualization, imagine that we only have 2-layers within the tissue. In other words, for each simulated photon, we have two partial paths, $L_1$ & $L_2$. Given the joint probability distribution $p(L_1, L_2)$, we would need to integrate over the line connecting $(L,0)$ to $(0,L)$ to obtain $I(L)$. 

ToF Filtering, on the other hand, simply refers to integrating $I(L)$ over some limit $L_{min}$ to $L_{max}$. Note that, setting the limits $L=0$ to $L=\inf$ emulates taking a CW intensity measurement. 

## Expressing Sensitivity in terms of ToF
We can rewrite the sensitivity equation to incorporate ToF related terms.
### Simplifying the Numerator 
From the sensitivity derivations,
```math
\frac {\delta I} {\delta TMP_{fetal}} = - \frac{\delta \mu_a} {\delta TMP_{fetal}} \times \sum I_iL_{fetal,i} \propto - \frac{\delta \mu_a} {\delta TMP_{fetal}} \times \frac{1}{N}\sum I_iL_{fetal,i}
```
Where, the summation is over every simulated photon, $i = 1, 2, ... N$.(N is in the range of 1e8 to 1e9 for our simulations)

The first term on the right-hand side is essentially a cosntant wrt to ToF filtering. So, we can focus on the second term for now.  

The second term can be rewritten as the expected value, $\mathbb{E}(I\times L_{fetal})$.  
```math
\mathbb{E}(I\times L_{fetal}) = \int I\times L_{fetal} p(L_{fetal}, I)dL_{fetal}dI
```
Note that the $I$ term here is the same as $I_{ToF}$ integrated over the entire range of $L$($L=0$ to $\infty$). Replacing the variables,
```math
\mathbb{E}(I\times L_{fetal}) =  \int_0^{\infty} I_{ToF}(L)L_{fetal}dL
```
Also not

### Replacing the $L_{fetal}$ term
The ToF measurements are always with respect to $L = \sum_j L_j$. It's very difficult to separate out the individual pathlengths from this term. However, under certain assumptions we can calculate a new term proportional to $L_{fetal}$ by using $I_{ToF}(L)$.  

Assume that we have the ToF measurements at two different points in time within a fetal pulsation. Note that this notion of time is not the same as the time axis in ToF measurements. Rather this corresponds to making two seprate sets of ToF measurements spaced out across time. We also assume these measurements are made such that the tissue layers do not change shape geometrically. We also assume that blood concentration within the maternal layer is the same between these measurements. Final assumption, the fetal layer pulsates slight allowing for a greater/lower blood concentration within the optical path, and a consequent change in the fetal layer $\mu_a$. Let's say due to this blood concentration change, the Hemoglobin concetration changes from $c_1$ to $c_2$. This leads to the absorption coefficient changing from $mu_{a,1}$ to $mu_{a,2}$.  

Taking the ToF measurements and dividing them point-by-point(at each L), we get
```math
R(L) = \frac{I_{1,ToF}(L)} {I_{2,ToF}(L)}
```
The difference between these two are on a single term, the $\mu_{a,fetal}$. Lets call them $\mu_{a,1}$ and $\mu_{a,2}$ for now. Then, $I_{2,ToF}(L)$ can be rewritten as 
```math
I_{2,ToF} = \int_{L_1} \int_{L_2} .. \int_{L_M} \mathfrak{1}_{(\sum_jL_j = L)} exp(-\sum_{j,j\neq fetal} \mu_{a,j}L_j)\frac{exp(-\mu_{a,2}L_{fetal})}{exp(-\mu_{a,1}L_{fetal})}exp(-\mu_{a,1}L_{fetal})p(L_1, L_2, .., L_M)
```
If we assume the average value $\hat L_{fetal}$, we can get a close approximate of $I_{2,ToF}$ as,
```math
I_{2,ToF} \approx I_{1,ToF} \times \frac{exp(-\mu_{a,2}\hat L_{fetal})}{exp(-\mu_{a,1}\hat L_{fetal})}
```
```math
R(L) = exp(-(\mu_{a,1} - \mu_{a,2})\hat L_{fetal}) = exp(-(c_1 - c_2)\epsilon_{Hb}(S)L_{fetal})
```
Where, $\epsilon_{Hb}(S)$ is the Hb extinction coefficient at the current fetal Saturation, S. If we choose the two measurement points close snough in time, this should be a constant. Since we assumed that there is no change in geometry due to pulsation, $c_1 - c_2$ is constant with respect to $L$. We can then simplify this relation as,
```math
|log(R(L))| \propto \hat L_{fetal}
```
Notice here that taking tha absolute value over $log(R(L))$ renders the **order** in which we chose the two points meaningless.  
(**NOTE:** This approximation might be invalid)

### Putting it All Together
We can write sensitivity as,
```math
Sensitivity \propto \frac{\int I_{ToF}(L) |log(R(L))| dL} {\int I_{ToF}(L) dL}
```
Where, the limits to both the integrals can be modified by ToF filtering. So after filtering, the sensitivity can be expressed as 
```math
Sensitivity_{filtered} = \frac{\int_{L_{min}}^{L_{max}} I(L) |log(R(L))| d(L)} {\int_{L_{min}}^{L_{max}} I(L) d(L)} 
```

## Optimizing Sensitivity
Our ToF filtering optimization goal now becomes finding the best limits for this integral. 
```math
argmax_{L_{min}, L_{max}} \frac{\int_{L_{min}}^{L_{max}} I(L) |log(R(L))| d(L)} {\int_{L_{min}}^{L_{max}} I(L) d(L)} 
```

## Optimization Comments
1. We want to maximize the numerator. The wider the integration range, the larger the integral becomes since all values are positive(As long as we consider there is no noise within the signal). In otherwords, we want to integrate over the widest range possible for the numerator.
2. We want to minimize the denominator. Again, as all terms are positive, the minimum value is 0. This is achieved when the upper and lower limit are equal. In otherwords, for the denominator, we want as narrow as a range as possible.
3. The optimum thus lies somewhere in between the narrowest and the widest range.
4. Since all the terms in the equation are avaialble to use from the ToF measurements, we can simply iterate over different combinations of integration ranges($L_1$ & $L_2$) and find the best range.
5. Non-Fetal Sensitive photons: We had a working definition this as photons whos $L_{fetal} = 0$. Notice the original sensitivity equation with the $L_{fetal}$ term in the numerator. Integrating over any region where $L_{fetal} = 0$ will not grow the numerator. However, it will grow the denominator, and consequently decrease sensitivity. Therefore, removing all Non-fetal sensitive photon would always show an improvement in sensitivity. How to filter these out photons out though is something I do not know. (One plausible idea might be to set some time threshold on the ToF measurement but even then this is not guaranteed to filter all such photons) 
