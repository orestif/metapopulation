# metapopulation
SIR model across multiple populations connected by migration

R source files (R Markdown files can be edited/compiled in RStudio).

For a complete description, see [Model Development](Model%20Development.html).

The three key source files that implement the three models (as of August 2015) are:

- "Metapop SIR adaptivetau.R": original model (based on Peel et al 2014) with Gaussian birth rate, constant death rate, frequency-dependent transmission (SIR), and explicit migration across a network of subpopulations.

- "Metapop DDD SIR adaptivetau.R": adds a second, density-dependent mortality term to prevent exponential increase in population sizes when birth pulses vary across the metapopulation.

- "Metapop DDD SIR adaptivetau.R": extension to two pathogen strains with no antigenic variations (perfect cross-protection).

All other files containing the word "series" have been used to run series of simulations from each of the three models, or analyse the resulting simulations. Note that I have not uploaded the output files from those simulation series as these are quite large.
