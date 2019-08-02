
%=========================================================================
%  Title : Ideal Adsorbed Solution Theory (IAST)
%  Author: Kaihang Shi
%          North Carolina State University, NC 27695, USA
%  Date  : 5/5/2019
%  
%  Model based on Myer & Prausnitz (1965)
%  This program will calculate the binary mixture adsorption isotherms 
%  over a pressure range based on the experimental pure isotherms
%=========================================================================

%=========================================================================
%  How to use:
%
%  1. Fit the pure adsorption/desorption isotherms using the proper model 
%  2. Modified the 'Input parameters' section below if necessary and enter
%  the fitted isotherm model parameters for each component
%  3. Press 'Run' button. Final results will be stored in the matrix 'sol'.
%
%  Note: working units are pressure [kPa], specific adsorbed amount [mmol/g]
%=========================================================================

% clean variable space
clear
clc

%=========================================================================
% Input parameters:
%
%  P_hi          upper bound of OUTPUT pressure limit in units of [kPa]
%  P_lo          lower bound of OUTPUT pressure limit in units of [kPa]
%  P_int         OUTPUT pressure interval in units of [kPa]
%  yi            mole fraction in gas phase
%  model_type    choose isotherm mode: 1-Langmuir, 2-Sips, 3-Toth
%  Lang_1        Langmuir isotherm parameter for species 1 [n_s,b]
%  Lang_2        Langmuir isotherm parameter for species 2
%  Sips_1        Sips isotherm parameter for species 1     [n_s,b,n]
%  Sips_2        Sips isotherm parameter for species 2
%  Toth_1        Toth isotherm parameter for species 1     [n_s,b,t]
%  Toth_2        Toth isotherm parameter for species 2
%=========================================================================
P_hi       = 100;         
P_lo       = 1;            
P_int      = 2;            
yi         = [0.15,0.85];    
model_type = 2;            
lang_1     = [0;0];        
lang_2     = [0;0];       
sips_1     = [3.801;7.309e-05;1.59];     
sips_2     = [0.8652;0.0009088;0.7573];     
toth_1     = [0;0;0];      
toth_2     = [0;0;0]; 

% Initialize parameters
sol = zeros(ceil((P_hi-P_lo)/P_int),7);
i = 0;

% Loop pressure
for press = P_lo:P_int:P_hi
    
    % Count
    i = i + 1;
    
    % Calculate partial pressure in the gas phase
    P1 = press*yi(1);
    P2 = press*yi(2);
    
    % Pick model
    % Langmuir model
    if model_type == 1
        fun1 = @(t)lang_1(1)*(lang_1(2)*t)./...
            (1+lang_1(2)*t)./t;
        fun2 = @(t)lang_2(1)*(lang_2(2)*t)./...
            (1+lang_2(2)*t)./t;
    % Sips model 
    elseif model_type == 2
        fun1 = @(t)sips_1(1)*(sips_1(2)*t).^(1/sips_1(3))./...
            (1+(sips_1(2)*t).^(1/sips_1(3)))./t;
        fun2 = @(t)sips_2(1)*(sips_2(2)*t).^(1/sips_2(3))./...
            (1+(sips_2(2)*t).^(1/sips_2(3)))./t;
    % Toth model
    elseif model_type ==3 
        fun1 = @(t)toth_1(1)*(toth_1(2)*t)./...
            (1+(toth_1(2)*t).^toth_1(3)).^(1/toth_1(3))./t;
        fun2 = @(t)toth_2(1)*(toth_2(2)*t)./...
            (1+(toth_2(2)*t).^toth_2(3)).^(1/toth_2(3))./t;
    else 
        error('Invalid model_type parameter!');
        
    end
    
    % 
    intfun = @(x)integral(fun1,0,P1/x)-integral(fun2,0,P2/(1-x));
    
    % Solve for adsorbed mole fraction x1, x2
    x1 = fzero(intfun,[0 1]);
    x2 = 1-x1;
    
    % P0 for each component where spreading pressure is equal
    P10 = P1/x1;
    P20 = P2/x2;
    
    % Calculate adsorption amount ni at P0
    % Langmuir model
    if model_type == 1
        n10 = lang_1(1)*(lang_1(2)*P10)/...
            (1+lang_1(2)*P10);
        n20 = lang_2(1)*(lang_2(2)*P20)/...
            (1+lang_2(2)*P20);
    % Sips model 
    elseif model_type == 2
        n10 = sips_1(1)*(sips_1(2)*P10)^(1/sips_1(3))/...
            (1+(sips_1(2)*P10)^(1/sips_1(3)));
        n20 = sips_2(1)*(sips_2(2)*P20)^(1/sips_2(3))/...
            (1+(sips_2(2)*P20)^(1/sips_2(3)));
    % Toth model
    elseif model_type ==3
        n10 = toth_1(1)*(toth_1(2)*P10)/...
            (1+(toth_1(2)*P10)^toth_1(3))^(1/toth_1(3));
        n20 = toth_2(1)*(toth_2(2)*P20)/...
            (1+(toth_2(2)*P20)^toth_2(3))^(1/toth_2(3));
    else
        error('Invalid model_type parameter!');
    end
    
    % Calculate total adsorbed mixture amount 
    ntot = 1/(x1/n10+x2/n20);
    % Adsorbed amount for component 1
    n1 = ntot*x1;
    % Adsorbed amount for component 2
    n2 = ntot*x2;
    
    % Final solution
    sol(i,1) = press;
    sol(i,2) = n1;
    sol(i,3) = n2;
    sol(i,4) = x1;
    sol(i,5) = x2;
    sol(i,6) = P10;
    sol(i,7) = P20;
    

end
