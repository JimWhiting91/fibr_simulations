# FIBR Simulations
Simulations of allele frequency dynamics among the FIBR guppy introduction experiments

### Author Contact Information
  * b.fraser@exeter.ac.uk
  * james.whiting@ucalgary.ca

### Usage and license information
If you use or are inspired by code in this repository please cite the following work or contact me about how to cite:

van der Zee et al. (20XX) *in review* [doi_link]()

Please also see the associated LICENSE.

---
### Environment Setup
All simulations were performed with [SLiM](https://messerlab.org/slim/) (v3.4), and analyses were performed with R (v4.0.3).

Clone this repo somewhere locally, and run the following to build the expected directory structure:
```
git clone git@github.com:JimWhiting91/fibr_simulations.git
cd fibr_simulations
mkdir outputs tables
```

### Simulation Workflow
  * The burn-in population is produced by running a population of 20,000 individuals until it reaches coalescence. This is done with the slim script `slim/guanapo_burnin.slim`.
  * Each FIBR introduction population is then run through its own demography based on census estimates. Demography files are available in `data/*_simulation_demography.txt`. These use the slim script `slim/simulate_fibr_introduction.slim`, and run using the bash wrapper `scripts/Run_fibr_AF_sims_interactive.sh`, which runs through each introduction population separately and loops through 200 iterations.
  * All analyses and figures are plotted using `R/fibr_sim_analysis_clean.R`

  * Additional genotype plotting is done with `R/fibr_genotype_plotting.R`
  * Per-generation census size estimates, based on the mean across ~8 months, are calculated using `R/00_make_demography_data_from_census.R`.

![Neutral AF Dynamics](./figs/FigureSX_neutral_simulation_results.png?raw=true "Neutral AF Dynamics")
