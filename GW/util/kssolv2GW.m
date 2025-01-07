function [GWinfor, optionsGW] = kssolv2GW(mol, options_in)
%   [GWinfor, optionsGW] = kssolv2GW(mol, options_in)
%   kssolv2GW: Converts the outputs of kssolv-scf to inputs required for
%              GW calculations, and sets up the parameters for the GW.
%
%   Inputs:
%       mol:       A @Molecule object from kssolv, containing molecular information.
%       options_in:A struct used to generate both the input data for the GW method
%                  and parameters used during GW calculations.
%                  Refer to the documentation for detailed descriptions.
%
%   Outputs:
%       GWinfor:   A @GWinfo object containing information required for GW calculations.
%                  Refer to the documentation for detailed descriptions.
%       optionsGW: A @GWOption object containing parameters for GW calculations.
%                  Refer to the documentation for detailed descriptions.
%
GWinfor = GWinfo();
optionsGW = GWOptions(mol, options_in);
GWinfor = kssolv2GW_info(GWinfor, mol, optionsGW.Groundstate);

end % end of function