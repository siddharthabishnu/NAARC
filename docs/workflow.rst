.. role:: bash(code)
   :language: bash

Workflow
========
This is an evolving document detailing the workflow to setup and run the NAARC configuration. 
At present it is very ad-hoc, but it's hoped this will be refined into a semi-automated
process.

LBCs
====

We are currently making use of the `bdy_msk` (at present undocumented in the NEMO manual) to define
the NAARC domain as a sub region of eORCA12. This enables us to avoid having to recompute the SBC 
weights files and ICs each time the domain boundaries are redefined.

- .ipynb to define :bash:`bdy_msk`, covering the Atlantic and Arctic reions, from :bash:`top_level` in :bash:`domain_cfg.nc`

  - :bash:`bdy_msk` is defined by 1s in the active domain and 0s elsewhere

- .ipynb to define :bash:`mask` from :bash:`bdy_msk` which is used as an input for pyBDY

  - :bash:`mask` is defined by 1s in the active domain, 0s on land and -1s elsewhere

- .ipynb to sub sample the :bash:`mask` in the region of each open boundary as pyBDY is not fully optimised for multiple LBCs over an extended domain.

  - in order to not confused pyBDY, the sub sampled :bash:`mask` must be modified by adding land points (0s) along the last row (see script for further details)



ICs
====

SBCs
====

DOM
===
