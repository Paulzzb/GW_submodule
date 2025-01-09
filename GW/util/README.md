This folder 
- transform inputs files into the inputs for GW module.
  - The code looks like [info_GW, options_GW] = ?2gw(...),
  -    where ? is the type of the input file (qe, ks, etc.).
  -    info_GW is the information of the system, class @GWInfo.
  -    options_GW is the options for GW module, class @GWOptions.

The converters do the following things:
1. Convert the output of the scf into ...

GW module accepts inputs in the format of @GWoptions and .
Currently, only inputs from kssolv are accepted.
# Structure of util folder
```
├── util
|   ├── kssolv2GW.m: Files are from kssolv ground-state calculation.
|   |    ├── kssolv2GW_info.m: Transform kssolv info file into GWinfo.
|   |    ├── kssolv2GW_initopt.m: Set GWOptions using information from kssolv.
|   |    ├── kssolv2GW_opt*.m: Subfunctions for setting GWOptions.
|   ├── qe2GW.m: Files are from Quantum Espresso ground-state calculation.
|                Currently, transform it into kssolv format, then to GWinput.
|                Not write yet...
|   ├── CoulombG: Calculate Coulomb in G-space.
|                 Need to be modified.
|   └── ...


```