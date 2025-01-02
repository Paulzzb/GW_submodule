### How to use the it.



## Files:
- demo_scf.m: run the scf iteration from KSSOLV 3.0.
- demo_GW.m: run the GW calculation.
- LiH222.m: Input file for the structure of moleculer.
  - Output: a @Moleculer class variable "mol".
- setGW.m: set GW parameters.
  - Output: a @GWOption class variable options_GW.
- setscf.m: set scf parameters.
  - Output: a "struct" variable "options_scf".

# Output files  
- GWoutput.mat: Energies obtained from GW calculation 
- molinfo.mat: a @Moleculer class variable "mol".
- scinfo.mat: result of scf iteration, which contains
  - mol: moleculer
  - H: Hamitonian
  - X0: wavefunctions
  - info: ...


#