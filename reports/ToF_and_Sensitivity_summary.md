# Relating Sensitivity to ToF Filters
## Defining Sensitivity 
We define sensitivity as the relative partial derivative of the intensity signal wrt any of the fetal TMP we care about(So either saturation or [Hb])
$$
Sensitivity = \frac {\delta I} {\delta TMP_{fetal}} / I 
$$

## ToF Measurements & ToF Filtering
A regular ToF measurement would give us an intensity at time T, $I(T)$. If we assume that the speed of light throughout the different layers within the tissue are equal, we can replace T with the total photon pathlength, L. This allows us to reparameterize $I(T)$ and break it down using Beer-Lambert's law
$$
I(T) = I(L) = exp(\sum_j \mu_{a,j}L_j)p(L) 
$$
Where, j denotes each tissue layer in the model, and $mu_{a,j}$ is its corresponding absorption cofficeint. $L_j$ denotes each photon's pathlength through only the j-th layer. $p(L) = p(L_1, L_2, ....L_{fetal})$ is the joint probability distribution over all the $L_j$. Also note that $L = \sum_j L_j$. 

We can visualize $I(L)$ as something akin to an $L_1$ norm wrt $L$. As an example, lets imagine there are only 2 layers. Then for each L, we can imagine a line connecting (0, L) $ (L, 0). Then the $p(L)$ in this case would include each point in the probability distribution on that line. (Such that $L_1 + L_2 = L$) 

On the other hand, ToF Filtering refers to integrating the ToF intensity($I(L)$) over some arbitrary range $L_1$ to $L_2$. Note that, setting $L_1 = 0$ and $L_2 = \inf$ emulates taking a CW measurement. 


## Simplifying the Numerator 
From the sensitivity derivations,
$$
\frac {\delta I} {\delta TMP_{fetal}} = - \frac{\delta \mu_a} {\delta TMP_{fetal}} \times \frac{1}{N}\sum I_iL_{fetal,i}
$$
Where, the summation is over every simulated photon, $i = 1, 2, ... N$.(N is in the range of 1e8 to 1e9 for our simulations)

The first term on the right-hand side is essentially a cosntant wrt to ToF filtering. So, we can focus on the second term for now.  

The second term is essentially the expected value, $\mathbb{E}(I\times L_{fetal})$.  
$$
\mathbb{E}(I\times L_{fetal}) = \int I\times L_{fetal} p(L_{fetal}, I)dL_{fetal}dI
$$
We can break down the $I$ term using Beer-Lambert's law.
$$
\int I\times L_{fetal} p(L_{fetal}, I)dL_{fetal}dI = \int exp(-\sum_j\mu_{a,j} L_j)L_{fetal} p(L)(dL) = \int I(L)L_{fetal}(dL)
$$
where, The equation is now integrated over all L. As in, $dL = dL_1 dL_2 ... dL_{fetal}$. __A ToF filter can effectively change the limits of this integral__.

## Replacing the $L_{fetal}$ term
The ToF measurements are always with respect to $L = sum_j L_j$. It's very difficult to separate out the individual pathlengths from this term. However, under certain assumptions we can calculate a new term proportional to $L_fetal$.  

Assume that we have the ToF measurements($I(L)$) at two different points in time within a fetal pulsation. Note that this notion of time is not the same as the time axis in ToF measurements. Rather this corresponds to making two seprate
sets of ToF measurements spaced out across time. We also assume these measurements are made such that the tissue layers do not change shape geometrically. We also assume that blood concentration within the maternal layer is the same between these measurements. Final assumption, the fetal layer pulsates slight allowing for a greater/lower blood concentration within the optical path, and a consequent change in the fetal layer $\mu_a$. Let's say due to this blood concentration change, the Hemoglobin concetration changes from $c_1$ to $c_2$. This leads to the absorption coefficient changing from $mu_{a,1}$ to $mu_{a,2}$.  
Taking the ToF measurements and dividing them point-by-point, we get
$$
R(L) = \frac{I_1(L)} {I_2(L)} = \frac{(exp(\sum_j \mu_{a,j}L_j)p(L))_1 }{(exp(\sum_j \mu_{a,j}L_j)p(L))_2}
$$ 
Now, realize that if the tissue geometry and the scattering coefficients do not change, $p(L)$ between these two points should also remain unchanged. Due the other assumptions we made, the only change will be on the fetal layer's $\mu_a$. Using these simplifications, the above equation can be rewriten as
$$
R(L) = exp(-(\mu_{a,1} - \mu_{a,2})L_{fetal}) = exp(-(c_1 - c_2)\epsilon_{Hb}(S)L_{fetal})
$$
Where, $\epsilon_{Hb}(S)$ is the Hb extinction coefficient at the current fetal Saturation, S. If we choose the two measurement points close snough in time, this should be a constant. With respect to the pathlengths, L, $c_1 - c_2$ is also a constant. We can then simplify this relation as,
$$
|log(R(L))| \propto L_{fetal}
$$
Notice here that taking tha absolute value over $log(R(L))$ renders the **order** in which we chose the two points meaningless.  

## Putting it All Together
We can write sensitivity as,
$$
Sensitivity \propto \frac{\int I(L) |log(R(L))| d(L)} {\int I(L) d(L)}
$$
Where, the limits to both the integrals can be modified by ToF filtering. So after filtering, the sensitivity can be expressed as 
$$
Sensitivity_{filtered} = \frac{\int_{L_1}^{L_2} I(L) |log(R(L))| d(L)} {\int{L_1}^{L_2} I(L) d(L)} 
$$

## Optimizing Sensitivity
Our ToF filtering optimization goal now becomes,
$$
L_1, L_2 = argmax_{L_1, L_2} \frac{\int_{L_1}^{L_2} I(L) |log(R(L))| d(L)} {\int{L_1}^{L_2} I(L) d(L)} 
$$

## Optimization Comments
1. We want to maximize the numerator. The wider the integration range, the larger the integral becomes since all values are positive(As long as we consider there is no noise within the signal). In otherwords, we want to integrate over the widest range possible for the numerator.
2. We want to minimize the denominator. Again, as all terms are positive, the minimum value is 0. This is achieved when the upper and lower limit are equal. In otherwords, for the denominator, we want as narrow as a range as possible.
3. The optimum thus lies somewhere in between the narrowest and the widest range.
4. Since all the terms in the equation are avaialble to use from the ToF measurements, we can simply iterate over different combinations of integration ranges($L_1$ & $L_2$) and find the best range.
5. Non-Fetal Sensitive photons: We had a working definition this as photons whos $L_{fetal} = 0$. Notice the original sensitivity equation with the $L_{fetal}$ term in the numerator. Integrating over any region where $L_{fetal} = 0$ will not grow the numerator. However, it will grow the denominator, and consequently decrease sensitivity. Therefore, removing all Non-fetal sensitive photon would always show an improvement in sensitivity. How to filter these out photons out though is something I do not know. (One plausible idea might be to set some time threshold on the ToF measurement but even then this is not guaranteed to filter all such photons) 
