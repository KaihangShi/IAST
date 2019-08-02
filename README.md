%=========================================================================
%  Title : Ideal Adsorbed Solution Theory (IAST)
%  Author: Kaihang Shi
%          North Carolina State University, NC 27695, USA
%  Date  : 5/5/2019
%  
%  Model based on Myer & Prausnitz (1965)
%  This program will predict the binary mixture adsorption isotherms 
%  based on the experimental pure isotherms
%=========================================================================

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
