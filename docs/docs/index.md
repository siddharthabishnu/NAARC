# Getting Started

**Welcome to the documentation for the North-Atlatic-ARCtic (NAARC) configuration**

## Introduction
As part of the Climate Change in the Arctic and North Atlantic Region and Impacts on the UK ([**CANARI**](https://canari.ac.uk)) project, the National Oceanography Centre has developed a regional NEMO configuration of the North Atlantic and Arctic oceans.

---

## Configuration :globe_with_meridians:

The NAARC configuration is based on NEMO v4.2.2 with a horizontal resolution of 1/12$^{\circ}$.

The key features are summarised below:

* 1/12$^{\circ}$ nominal horizontal resolution (j=3605, i=4320).
* 75 vertical Multi-Evelope s-levels (MEs).
* Coupled to [SI$^{3}$](https://doi.org/10.5281/zenodo.7534900) sea ice engine.
* Initialised from [eORCA025 NOC-Near-Present-Day](https://noc-msm.github.io/NOC_Near_Present_Day) (1979) simulations.
* Forced with JRA55-do (v1; 1976-2023) atmospheric and riverine forcing.

For more details on the NAARC configuration see [Deep Dives: Model Configurations].

[Deep Dives: Model Configurations]: deep_dives.md#model-configurations

---

## Quick Start :rocket:

### Installation

To get started, check out and set up an instance of the NAARC GitHub [repository](https://github.com/NOC-MSM/NAARC):

```sh
git clone git@github.com:NOC-MSM/NAARC.git
```

??? tip "Helpful Tip..."

    * **It is not advised to checkout the respository in your home directory.**

Next, run the setup script to download [NEMO](https://www.nemo-ocean.eu) & compile the tools and configurations. At present the setup 
script has been tested and will checkout, compile and run the NAARC (NEMO 4.2.2) code on: ARCHER2 for Cray-MPICH and GNU-MPICH, and Anemone for iFort.

=== "Anemone"
    ```sh
    cd NAARC

    ./NAARC/scripts/setup/NAARC_setup -p $PWD/NAARC_RUNS  -r $PWD/NAARC -n 4.2.2 -x 2 -m anemone -a impi -c ifort
    cd NAARC_RUNS/nemo/cfgs/NAARC//
    cp -rP EXPREF EXP_MYRUN
    cd EXP_MYRUN
    ln -s ../INPUTS/domain_cfg_mes.nc domain_cfg.nc
    ```

=== "Archer2"
    ```sh
    cd NAARC

    ./NAARC/scripts/setup/NAARC_setup -p $PWD/NAARC_RUNS  -r $PWD/NAARC -n 4.2.2 -x 2 -m archer2 -a mpich -c gnu
    cd NAARC_RUNS/nemo/cfgs/NAARC/
    cp -rP EXPREF EXP_MYRUN
    cd EXP_MYRUN
    ln -s ../INPUTS/domain_cfg_mes.nc domain_cfg.nc
    ```

In the case of `-m archer2` one can also complile with `-c cray`.

### Running An Experiment

Edit the project code and options in  `runscript.slurm` then:

```
sbatch runscript.slurm -y 1979 -s 1
```

This will produce a 5 day mean output from the beginning of 1979. The run should take 15 minutes to complete once in the machine.
To extend this to a full year's simulation, edit the `runscript.slurm` file and change:

```
   nn_itend    =    8928   !  last  time step (std 5840)
!  nn_itend    =    XXX_TEN_XXX   !  last  time step (std 5840)
```
to
```
!  nn_itend    =    8928   !  last  time step (std 5840)
   nn_itend    =    XXX_TEN_XXX   !  last  time step (std 5840)
```

### Forcing data:

[NAARC](https://gws-access.jasmin.ac.uk/public/jmmp/NAARC/)

_this is automatically transferred when the setup script is executed_

For ARCHER2 users these data are held under `/work/n01/shared/NAARC` and `/work/n01/shared/nemo/FORCING` and are linked during the setup.
