# Adaptation-to-TemporalChange
R scripts for simulations of evolution of adaptation and plasticity to temporal environmental change

Journal article: "Predicting the evolution of adaptation and plasticity from temporal environmental change" Gallegos et al. 2023

bioRxiv preprint: https://doi.org/10.1101/2023.02.12.528221

Shinny app: https://posit.cloud/content/6515271

DOI: 10.5281/zenodo.14029488

Contents:

- Evolution of environmental tolerance
 - Function to simulate environments (sim.env)
 - Function to simulate evolution under each environment (sim.evo)
 - Function to extract mean breeding value, mean plasticity, and tolerance curves at equilibrium (evo.out)
 - Function to plot tolerance curve at equilibrium with variation in evolved plasticity (breadth) (plot.tolcurve)
 - Simulate evolution under baseline scenario and plot outcomes
 - Simulate evolution under each environmental change scenario and plot outcomes

- Environmental predictability and evolution of plasticity (simulate evolution under each environment and extract and plot evolved mean plasticities at equilibrium)

- Case study
 - Download sea surface temperature data, crop to southeast Australia, derive temperature variables, and plot maps
 - Extract temperature time series for each location, simulate evolution of thermal tolerance for each location, and plot outcomes
