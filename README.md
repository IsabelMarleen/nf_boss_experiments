# nf_boss_experiments
This repository contains a Nextflow workflow to run different simulated experiments with BOSS-RUNS. The following workflow is implemented with more under development (see branch parameter_sweep):

1. **Human trio sequencing**: The trio sequencing workflow creates consensus genomes using a reference and vcf calls. From these consensus genomes, fragments are generated using [NanoSim](https://github.com/bcgsc/NanoSim.git) to simulate ONT like errors and then barcoded. These fragments undergo some preprocessing (reimplementation of [prepare_simulation_data.smk](https://github.com/goldman-gp-ebi/BOSS-RUNS/blob/main/scripts/prepare_simulation_data.smk)). They are then sequenced *in silico* and analysed.


## Repository structure
This repository is inspired by nf-core standards, though not fully consistent with their guidelines.
* The ```bin``` folder contains python and R scripts needed by the different modules.
* The ```conf``` folder contains configuration files for the pipeline (base.config), the modules (modules.config) and a test run (test.config).
* The ```envs``` folder contains environment specification files to build conda environments needed by the pipeline.
* The ```modules``` folder contains all the individual processes written for this pipeline.
* The ```workflows``` folder contains workflow files defining larger subprocesses within the pipeline.
* ```main.nf``` is the core nextflow file and defines the highest level workflow. This file is likely to change as the range of workflows possible is expanded.
* ```nextflow.config``` is the top-level configuration file and describes the possible parameters and other fundamental options.

## Requirements
* Nextflow
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

**NOTE** Currently, the testing profile requires a precomputed profile for simulating fragments. This will work automatically, if you clone the nanosim repository one directory level above this repository.

## Usage
### Setup
To run your own experiments, you will need to specify a number of input parameters. I recommend writing a parameter file for each of your experiments. An example can be found in ```conf/example-params.json```. Enabling singularity will improve the speed and performance as the initial building of conda environments can take a long time. Providing a container for all processes is a work in progress.

### Running the pipeline
The pipeline can be run by executing

```nextflow run main.nf [-profile ... -c ...] -with-conda -params-file custom_params.json```

