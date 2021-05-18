#!/bin/python

# Python script to simulate a neutral burn-in for GH using msprime
# GH is assumed to have an Ne of 20,000 from Whiting et al. 2021 PLOS Genetics
# Mutation rate is assumed at 4.89e-8 from Kunster et al. 2016 PLOS One

# assumes conda activate msprime-env

import msprime, pyslim

# Set variables
GH_Ne = 20000
chrom_size = 1e6
mu = 4.89e-8

# Do simulation
ts = msprime.simulate(sample_size=GH_Ne*2, Ne=GH_Ne, length=chrom_size,
    mutation_rate=mu, recombination_rate=1e-8)
slim_ts = pyslim.annotate_defaults(ts, model_type="WF", slim_generation=1)
slim_ts.dump("outputs/GH_burnin_msprime.trees")
