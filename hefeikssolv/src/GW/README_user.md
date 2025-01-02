# GW code

## Current progress

- [ ] GW_GPP
  - [x] Move and adjust from _swp to here
- [x] GW_FULLFREQ_CD
  - [x] Move and adjust from _swp to here
  - [x] Some adjustmwhos
  - [ ] ents are needed, check gw_fullfreq_cd for details
- [ ] GW_FULLFREQ_RA
  - [ ] Move and adjust from _swp to here
  - [ ] Some details work are needed, including numerical integral analysis.

## To do
- [] finish and test gw_fullfreq_cd_, check if it is the same as gw_fullfreq_cd.

- [] Remove redundant files.
- [] (For others to do) GPU support on GW module, mainly isdf_kernelg.m
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Structure of GW module

- $KS_HOME/src/GW
  - Common: contains commonly used function.
    - Cauchy integral
    - Fourier transform
    - orbital pair functions
    - find corresponding g index
  - @ksinfo (See 'classes in GW, @ksinfo' for details)
  - @gvec (See 'classes in GW, @gvec' for details)
  - options (See 'classes in GW, GWOptions' for details)
  - testGW/: testfunctions in GW module, for both users and developers.
  - isdf: A stable version for our GW code.
# Workflow of GW module

mol, options_in --(gwCalculation.m)--> GWoutput, which contains

1. mol, options_in --(GWOptions.m)--> options
   - This step contains initialization of inputs, and takes little time.
2. mol, options.Groundstate --(ksinfo.m)--> ksinfo
   - This step contains
     - calculating groundstate information.
     - Preparing some other datas for GW module.  
   - This step costs the same as normal groundstate calculation using scf.
3. ksinfo, options --(gw_*.m)--> GWoutput
   - This step is the main part of GW module.
   - Normally it takes O(n_e^3) times
     - ~ 50 seconds for Si64 with COHSEX if using k-means ISDF.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Samples

There are examples named $home/test/testGW/testGW_*.m, and they all follow the steps below.

- Step 1: Generate 'mol' in Moleculer class.
  - Examples of constructing 'mol' are in \$home/example.
- Step 2: Generate a struct named 'options_in' and set the field values.
  - Examples are in \$home/test/testGW/testGW_\*.m and
  \$home/test/testGW/testGW_\*.m.
  - Please refer to the section 'Classes in GW module/GWOptions'
    for the specific meanings of each field in options_in.
- Step 3: Call gwCalculation.m.
<!-- 
- Step 3: Use 'GWOptions' with 'mol' and 'options_in' as inputs to generate 
  'options'.
  - Please refer to the section 'Classes in GW module/GWOptions' for details.
- Step 4: Use class 'ksinfo' with 'mol' and 'options.Groundstate' as inputs
  to implement scf and generate 'ksinfor'.
  - Please refer to the section 'Classes in GW module/ksinfo' for details.
- Step 5: Use gwCalculation.m with 'ksinfor' and 'options' to do GW calculation. 
-->
# Classes in GW module

__WARNING: WE STRONGLY SUGGEST THE USERS NOT TO CHANGES VALUES OF FIELDS IN THE
FOLLOWING CLASSES DIRECTLY.__
In the following of this section, we use
__XX1 <-- XX2 (value1): info__
where

- __XX1__: field in Outputs
- __XX2__: field in inputs, which set the values of XX1
- __value1__: default value if XX2 is not the field of Inputs
- __info__: means of XX1.

1. GWOptions
A class to init and reschedule all parameters for GW calculation.
   - Inputs: options_in, check testGW/testGW_cohsex.m for examples.
   - Outputs: GWOptions class 'options', with the fields as followed
     - Groundstate: For groundstate calculation.
       - isGW <-- isGW (true): Do GW or not.
       - isBSE(Not support now) <-- isBSE (false): Do BSE or not.
       - frequency_dependence <-- frequencydependence (0)
       - nv <-- nv (mol.nel / 2 if spin degenerate): number of unoccupied states
       - nc <-- nc (mol.nel / 2 if spin degenerate): number of occupied states
       - amin <-- amin (5): Used when coulomb_truncation=2 (spherical)
       - coulomb_truncation <--  coulomb_truncation (2):
         - 2: spherical truncation
         - 6 (not tested yet): cell slab truncation.
       - options_kssolv <-- options_kssolv (none):
         - options for running scf, default setting in @ksinfo/gwsetup.m  
       - input <-- input ('kssolv'): where the inputs file from
       - inputfile <-- inputfile (none): read scf result if not empty.
     - Constant: Constants for GW calculation
       - __nv__ <-- nv (mol.ne / 2 if spin degenerate): number of all occupied
         states
       - __nc__ <-- nc (mol.ne / 2 if spin degenerate): number of all unoccupied
         states taken for Groundstate calculation
       - __nv_ener__ <-- nv_ener (mol.ne / 2 if spin degenerate):
         number of states under fermi level that we calculate GW energy.
       - __nc_ener__ <-- nc_ener (mol.ne / 2 if spin degenerate):
         number of states above fermi level that we calculate GW energy.
       - __nv_oper__ <-- nv_oper (mol.ne / 2 if spin degenerate):
         number of states under fermi level for constructing $\chi$.
       - __nc_oper__ <-- nc_oper (mol.ne / 2 if spin degenerate):
         number of states above fermi level for constructing $\chi$.
       - __DO NOT CHANGE VALUES IN CONSTANT__
     - ISDFCauchy
       - __isISDF__ <-- isISDF (false): Apply ISDF to calculation or not.
         - we do not support gw_gpp with isISDF == true.
       - __isCauchy__ <-- isCauchy (false): Calculate C*Omega^{-1}*C'
         with Cauchy integral or not
       - optionsCauchy <-- optionsCauchy
         - froErr <-- optionsCauchy.froErr (1e-6): integral tolerate.
         - MaxIter <-- optionsCauchy.froErr (10): number of the Max iteration
           for numerical integral.
       - __vcrank_ratio__ <-- vcrank_ratio (8): ratio for decomposite M_{vc}.
       - __vsrank_ratio__ <-- vsrank_ratio (8): ratio for decomposite M_{vs}.
       - __ssrank_ratio__ <-- ssrank_ratio (8): ratio for decomposite M_{ss}.
       - __The following options for ISDF goes ISDF.m for details.__
       - exxmethod <-- optionsISDF.exxmethod ('qrcp')
       - (Some other unimportant ISDF parameters, go ISDF module for details.)
     - GWCal
       - __fileName__ <-- fileName ('GW_output.mat'): file to save GW energy.
       - freq_dep <-- frequency_dependence (0):
         - 0: COHSEX
         - 1: GPP
         - 2: FULL FREQUENCY.
       <!-- - <-- (): -->
       <!-- - <-- (): -->
       - ... (Others are for gpp and full-frequency)
2. ksinfo
   - Input: mol, options.Groundstate
   - Output: Groundstate informations for further calculation.
     - mol:
     - coulG: size(ng) * 4, each row contain [g_1, g_2, g_3, v([g])]
       - The coulG contains the G-grid in values related with 'ksinfo.gvec', not the G-grid in "ksinfo.gvec2" and "ksinfo.gvecrho"
     - coulG0: correction on g = 0.
     - bdot:
     - qk: q points, not used yet.
     - ntot: number of discrete points in real space.
     - vol: volume of the system
     - ne: number of electron.
     - nv: number of occupied states
     - gvec, gvec2, gvecrho: @gvec class variable.
       - gvec: for main GW calculation.
       - gvec2: not used yet.
       - gvecrho: gvec related to rho (since in gpp model we takes larger ecut for rho.)
     - Vxc: exchange and corrolation energies for states.
     - idxnz: equal to gvec.idxnz (We take it out since it is commonly used in GW module.)
     - aqs: precalculated orbital pair functions.
       - Defaultly, we do not precalculate and save it.
   - __DO NOT CHANGE VALUES IN THIS CLASS__
3. gvec
   save the grid info in fourier space.
   - Input: mol, varargin
     - __ecut__ (5): sphere truncation energy in reciprocal space,
       __Ha as unit__
     - method ('kssolv'): where the input files from, only support 'kssolv' now.
     - dir (None): if choose other method, where the input files from.
