# nf_boss_experiments
This repository contains a Nextflow workflow to run different simulated experiments with BOSS-RUNS. It combines the creation of simulated reads based on benchmark vcf calls using Nanosim, the preprocessing needed to run a BOSS-RUNS simulated run (reimplementation of [prepare_simulation_data.smk](https://github.com/goldman-gp-ebi/BOSS-RUNS/blob/main/scripts/prepare_simulation_data.smk)) the sequencing and then the analysis of the run (reimplementation of the [BOSS_Simulation_Example](https://github.com/W-L/boss_simulation_example.git) repository.

## Setup
This pipeline requires the [BOSS-RUNS](https://github.com/goldman-gp-ebi/BOSS-RUNS.git) and [NanoSim](https://github.com/bcgsc/NanoSim.git) repositories to be installed. It also requires a number of other programmes ideally in separate environments as managed by conda, documentation on this to follow.
