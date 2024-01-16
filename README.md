# TFO Sensitivity


How sensitive are our simulated TFO signals with respect to each of the tissue model parameters.

 We determine the partial derivative of our output signal w.r.t. these TMPs using both numerical and analytical formulae. The derivations for the analytical parts are in reports. We also care about the comparative sensitivity between different modelling schemes.

## Project Organization
    ├── LICENSE
    ├── README.md          <- The top-level README for developers using this project.
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    ├── notebooks          <- 
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── tfo_sensitivity    <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── calculate_intensity  <- Intensity calculation related scripts
    │   │
    │   ├── jacobian       <- Calculating Jacobians (both numerically & analytically)
    │   │
    │   ├── data           <- loading/saving data
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations

## Jacobian



