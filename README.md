# nf_boss_experiments
This repository contains a Nextflow workflow to run different simulated experiments with BOSS-RUNS. It combines the creation of simulated reads based on benchmark vcf calls using Nanosim, the preprocessing needed to run a BOSS-RUNS simulated run (reimplementation of [prepare_simulation_data.smk](https://github.com/goldman-gp-ebi/BOSS-RUNS/blob/main/scripts/prepare_simulation_data.smk)) the sequencing and then the analysis of the run (reimplementation of the [BOSS_Simulation_Example](https://github.com/W-L/boss_simulation_example.git) repository. Each of the four subprocesses (creation of input genome, preprocessing, sequence simulation and analysis) is handled by a separate nextflow module and can be executed individually.

## Repository structure

## Requirements
* Nextflow
* [BOSS-RUNS](https://github.com/goldman-gp-ebi/BOSS-RUNS.git)
* [NanoSim](https://github.com/bcgsc/NanoSim.git)

## Installation
Clone this repository using

`git clone https://github.com/IsabelMarleen/nf_boss_experiments.git`

Then follow the steps below under Usage.

## Usage
### Setup
Running this workflow requires
### Running the pipeline


