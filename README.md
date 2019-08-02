This is a MATLAB code for Ideal Adsorbed Solution Theory (IAST) calculations. 
The IAST applied in this code is based on the original model by Myer & Prausnitz (1965)

The IAST is a theory to predict the binary mixture adsorption isotherms based on the 
pure adsorption isotherm data.


%=========================================================================
%  How to use:
%
%  1. Fit the pure adsorption/desorption isotherms using the proper model 
%  2. Modified the 'Input parameters' section if necessary and enter
%     the fitted isotherm model parameters for each component
%  3. Press 'Run' button in MATLAB. Final results will be stored in the matrix 'sol'.
%
%  Note: working units are pressure [kPa], specific adsorbed amount [mmol/g]
%=========================================================================
