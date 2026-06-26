# nf_boss_experiments
This repository contains a Nextflow workflow to run different end-to-end simulated sequencing experiments with [BOSS-RUNS](https://github.com/goldman-gp-ebi/BOSS-RUNS.git).

Depending on the type of experiment chosen, a combination of up to five subworkflows is executed. These subworkflows perform the following tasks in this order:
* VARIANT_INPUT: If no individual genome is available, this subworkflow can create it using `bcftools consensus` from variant calls (vcf files) and the reference genome they were called on.
* SIMULATE_FRAGMENTS_BOSS: The individual genome(s) created in the first step or provided directly by the user are used to create simulated reads with ONT-like error profiles using [NanoSim](https://github.com/bcgsc/NanoSim.git).
* PREPROCESS_BOSS: The reads created in the previous step or provided directly by the user undergo preprocessing required by the sequencing simulation, such as the creation of truncated reads, mapping of reads to the reference and creating the necessary config files for the sequencing. This step started as a reimplementation of [prepare_simulation_data.smk](https://github.com/goldman-gp-ebi/BOSS-RUNS/blob/main/scripts/prepare_simulation_data.smk) but has since diverged.
* SEQUENCE_PROFILE_BOSS: This subworkflow uses the files created in the previous step and uses the internal simulation mode in the BOSS-RUNS software to sequence all available reads for two experiment regime, one using dynamic adaptive sampling (boss) and the second a control. This subworkflow also creates time and memory profiles of the sequencing run but alternatives without profiling can be made available.
* ANALYSE_BOSS: The output files and logs produced by the sequencing simulation and analyses the run in regards to how for example rejections, coverage and evenness change over sequencing time. This step started as a reimplementation of the [BOSS_Simulation_Example](https://github.com/W-L/boss_simulation_example.git) repository but has since diverged.

## Description of available workflows
In progress

## Repository structure
The repository structure is oriented around nf-core guidelines and looks as follows:
* `bin/`: the binary directory contains scripts needed by processes and is automatically added to the PATH during execution
* `conf/`: the configuration directory contains a number of configuration files for the basic setup (base.config), module setup (modules.config), the test profile (test.config) and and an example parameter json file (example-params.json)
* `envs/`: The environment directory contains yaml files for setting up conda environments needed by some processes.
* `modules/`: The modules directory contains a subdirectory for all local module files (`modules/local/`) and one for imported nf-core modules (`modules/nf-core/`). Module files are nextflow files that contain the specification of a single process.
* `subworkflows/`: The subworkflow directory contains local and nf-core imported subworkflows in their respective subdirectories. A subworkflow consists of a series of processes. It takes one or more input channels and produces one or more output channels.
* `workflows/`: The workflow directory contains local and nf-core imported workflows in their respective subdirectories. A workflow consists of a series of subworkflows and processes. Unlike subworkflows, workflows create their inputs from input parameters and can thus be executed independently.
* `main.nf`: This is the main nextflow file that controls code flow and calls different named workflows.
* `nextflow_schema.json`: This file is a work in progress, where all potential input parameters and their types are defined. It is necessary for implementing parameter validation, which is under development.
* `nextflow.config`: This config file defines input parameters, the pipeline manifest and other configuration options.

## Requirements
* Nextflow==25.04.6
* Singularity
* Conda

## Installation
Clone this repository using

`git clone https://github.com/IsabelMarleen/nf_boss_experiments.git`

Then follow the steps below under Usage.

## Usage
### Checking the installation
To check whether the installation was successful you can run the test setup. This feature is still experimental and not quite as minimal as desired. The full execution takes about ten minutes. It is run as follows:

`nextflow run main.nf -profile test,<additional profiles> -c <custom_config_files>`

### Setting up an experiment
Given the large number of parameters, I recommend creating a json file with the parameters (see `conf/example-params.json` for an example). The implementation of parameter validation is ongoing.

Additionally, you might want to create a custom nextflow profile to set details about your specific compute environment (see [nextflow documentation](https://docs.seqera.io/nextflow/config#config-profiles) for details) or a custom configuration file to configure custom resource allocations for individual processes beyond the defaults set in `conf/modules.config` (see [nextflow documentation](https://docs.seqera.io/nextflow/config) for details).

Make sure that singularity and conda are enabled in your compute environment. Extending the pipeline to docker is in progress, but not yet implemented. Similarly, another goal is to make the pipeline runnable with just singularity or conda.

### Running the pipeline
After setting up an experiment, the pipeline can be run as follows:

`nextflow run <path to repository>/main.nf -c <custom.config> -c <manual_resources.config> -with-conda -params-file <params.json>`