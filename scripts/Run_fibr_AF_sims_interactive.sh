#!/bin/bash
# Run in parallel

# Set environment
MAIN=~/Exeter/fibr_simulations
cd $MAIN

# Set up metadata
fibr_pops=(CA TY LL UL)
fibr_start_sizes=(104 128 76 76)

# Make somewhere to store results
NO_SELECTION=outputs/FIBR_runs_without_selection
mkdir $NO_SELECTION

for i in {1..200}
do
    echo "STARTING ITER $i"

    slim -d founding_size="${fibr_start_sizes[0]}" -d 'demo_data="data/CA_simulation_demography.txt"' -d 'burnin_path="data/GH_burnin.txt"' slim/simulate_fibr_introduction.slim > $NO_SELECTION/CA_test_run_iter${i}.res

    slim -d founding_size="${fibr_start_sizes[1]}" -d 'demo_data="data/TY_simulation_demography.txt"' -d 'burnin_path="data/GH_burnin.txt"' slim/simulate_fibr_introduction.slim > $NO_SELECTION/TY_test_run_iter${i}.res

    slim -d founding_size="${fibr_start_sizes[2]}" -d 'demo_data="data/LL_simulation_demography.txt"' -d 'burnin_path="data/GH_burnin.txt"' slim/simulate_fibr_introduction.slim > $NO_SELECTION/LL_test_run_iter${i}.res

    slim -d founding_size="${fibr_start_sizes[3]}" -d 'demo_data="data/UL_simulation_demography.txt"' -d 'burnin_path="data/GH_burnin.txt"' slim/simulate_fibr_introduction.slim > $NO_SELECTION/UL_test_run_iter${i}.res
done
