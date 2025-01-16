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

 - in order to not confuse pyBDY, the sub sampled :bash:`mask` must be modified by adding land points (0s) along the last row (see script for further details)

- once the mask files are defined pyBDY is used to generate the :code:`coordinates.bdy.nc` files for the Atlantic and Pacific boundaries

 - for our reanalysis simulations the NOC Near Present Day (NPD) JRA forced eORCA025 simluation is used a source data for LBCs
 - pyBDY mesh-mask information 

- as the ice/snow outputs from the NPD run are not in the correct form for NEMO LBC inputs we have to divide :code:`sivolu` and :code:`snvolu` by :code:`siconc` to get the required :code:`sithic` and :code:`snthic`

- due to the nature on the interpolation between the source NPD data and the LBCs there may also be a difference between the total volume transport across the boundaries when comparing the source data and boundary files we are generating. pyBDY currectly doesn't allow for any transport adjustment, so this is done in an .ipynb as a post processing step.

- the other thing that may happen when generating the T/S boundary conditions is that the water column becomes statically unstable. Again, we handle this as a postprocessing step mixing down any denser water an unstable water column.


ICs
====

SBCs
====

DOM
===
