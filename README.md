# nf_boss_experiments
This repository contains a Nextflow workflow to run different end-to-end simulated sequencing experiments with [BOSS-RUNS](https://github.com/goldman-gp-ebi/BOSS-RUNS.git).

Depending on the type of experiment chosen, a combination of up to five subworkflows is executed. These subworkflows perform the following tasks in this order:
* **VARIANT_INPUT**: If no individual genome is available, this subworkflow can create it using `bcftools consensus` from variant calls (vcf files) and the reference genome they were called on.
* **SIMULATE_FRAGMENTS_BOSS**: The individual genome(s) created in the first step or provided directly by the user are used to create simulated reads with ONT-like error profiles using [NanoSim](https://github.com/bcgsc/NanoSim.git).
* **PREPROCESS_BOSS**: The reads created in the previous step or provided directly by the user undergo preprocessing required by the sequencing simulation, such as the creation of truncated reads, mapping of reads to the reference and creating the necessary config files for the sequencing. This step started as a reimplementation of [prepare_simulation_data.smk](https://github.com/goldman-gp-ebi/BOSS-RUNS/blob/main/scripts/prepare_simulation_data.smk) but has since diverged.
* **SEQUENCE_PROFILE_BOSS**: This subworkflow uses the files created in the previous step and uses the internal simulation mode in the BOSS-RUNS software to sequence all available reads for two experiment regime, one using dynamic adaptive sampling (boss) and the second a control. This subworkflow also creates time and memory profiles of the sequencing run but alternatives without profiling can be made available.
* **ANALYSE_BOSS**: The output files and logs produced by the sequencing simulation and analyses the run in regards to how for example rejections, coverage and evenness change over sequencing time. This step started as a reimplementation of the [BOSS_Simulation_Example](https://github.com/W-L/boss_simulation_example.git) repository but has since diverged.

## Description of available workflows
The following workflow is implemented with more under development:

* **Human trio sequencing**: The trio sequencing workflow creates consensus genomes using a reference and vcf calls. From these consensus genomes, fragments are generated using NanoSim to simulate ONT like errors and, if specified in the input parameters, barcoded. These fragments undergo preprocessing, then *in silico* sequencing and analysis.

## Repository structure
The repository structure is oriented around nf-core guidelines and looks as follows:
* The ```bin/``` directory contains python and R scripts needed by the different modules and is automatically added to the PATH during execution
* The ```conf/``` directory contains configuration files for the pipeline (base.config), the modules (modules.config) and a test run (test.config).Furthermore, it contains an example parameter json file (example-params.json).
* The ```envs/``` directory contains environment specification files to build conda environments needed by the pipeline.
* The ```modules/``` directory contains a subdirectory for all local module files (`modules/local/`) and one for imported nf-core modules (`modules/nf-core/`). Module files are nextflow files that contain the specification of a single process.
* `subworkflows/`: The subworkflow directory contains local and nf-core imported subworkflows in their respective subdirectories. A subworkflow consists of a series of processes. It takes one or more input channels and produces one or more output channels.
* The ```workflows/``` directory contains local and nf-core imported workflows in their respective subdirectories. A workflow consists of a series of subworkflows and processes. Unlike subworkflows, workflows create their inputs from input parameters and can thus be executed independently.
* `main.nf` is the core nextflow file that controls code flow and calls different named workflows.
* `nextflow_schema.json` is file, where all potential input parameters and their types are defined. It is necessary for implementing parameter validation which is under development.
* `nextflow.config` is the top-level configuration file and describes the possible parameters and other fundamental options.


## Requirements
* Nextflow==25.04.6
* Conda
* Most processes also provide singularity containers instead of conda

## Installation
Clone this repository using

`git clone https://github.com/IsabelMarleen/nf_boss_experiments.git`

Then follow the steps below under Usage.

### Testing
Once you have installed the repository and the other requirements, you can test the pipeline on the commandline:
```nextflow run main.nf -profile test```.

If you have an institutional profile or any additional custom configurations like where to write the work directory or the output to, then you can add additional profile and configuration files using ```-profile``` and the ```-c``` flags respectively.

The test run takes about ten minutes in an hpc environment.

**NOTE** Currently, the testing profile requires a precomputed profile for simulating fragments with NanoSim. This will work automatically, if you clone the nanosim repository one directory level above this repository.

## Usage
### Setting up an experiment
Given the large number of parameters, I recommend creating a json file with the parameters (see `conf/example-params.json` for an example). 
To run your own experiments, you will need to specify a number of input parameters. I recommend writing a parameter file for each of your experiments. An example can be found in ```conf/example-params.json```. The implementation of parameter validation is ongoing.

Additionally, you might want to create a custom nextflow profile to set details about your specific compute environment (see [nextflow documentation](https://docs.seqera.io/nextflow/config#config-profiles) for details) or a custom configuration file to configure custom resource allocations for individual processes beyond the defaults set in `conf/modules.config` (see [nextflow documentation](https://docs.seqera.io/nextflow/config) for details). Enabling singularity will improve the speed and performance as the initial building of conda environments can take a long time. Providing a singularity and docker container for all processes is a work in progress.

### Running the pipeline
After setting up an experiment, the pipeline can be run as follows:

`nextflow run <path to repository>/main.nf [-c <custom.config> -c <manual_resources.config>] -with-conda -params-file <params.json>`