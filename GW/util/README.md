This folder 
- transform inputs files into the inputs for GW module.
  - The code looks like [info_GW, options_GW] = ?2gw(...),
  -    where ? is the type of the input file (qe, ks, etc.).
  -    info_GW is the information of the system, class @GWInfo.
  -    options_GW is the options for GW module, class @GWOptions.

GW module accepts inputs in the format of @GWoptions and .
Currently, only inputs from kssolv are accepted.
# Structure of util folder
```
├── util
|   ├── kssolv2GW.m: Files are from kssolv ground-state calculation.
|   ├── qe2GW.m: Files are from Quantum Espresso ground-state calculation.
|                Not write yet...
|   └── ...


```