% Digital Earth Lab
% www.DigitalEarthLab.com
% Written by Andrew Pethick 2013
% Last Updated October 29th 2013
% Licensed under WTFPL

function [apparentResistivity, phase] = modelMT(resistivities, thicknesses,period)

for k = 1 : length(period)
mu = 4*pi*1E-7; %Magnetic Permeability (H/m)  
w = 2 * pi / period(k); %Angular Frequency (Radians);
n=length(resistivities); %Number of Layers

impedances = zeros(n,1);
%Layering in this format
%  Layer     j
% Layer 1    1
% Layer 2    2
% Layer 3    3
% Layer 4    4
% Basement   5

% Steps for modelling (for each geoelectric model and frequency)
% 1. Compute basement impedance Zn using sqrt((w * mu * resistivity))
% 2. Iterate from bottom layer to top(not the basement)
    % 2.1. Calculate induction parameters
    % 2.2. Calculate Exponential factor from intrinsic impedance
    % 2.3 Calculate reflection coeficient using current layer
    %          intrinsic impedance and the below layer impedance
    
% 3. Compute apparent resistivity from top layer impedance
        %   apparent resistivity = (Zn^2)/(mu * w)

%Symbols
% Zn - Basement Impedance
% Zi - Layer Impedance
% wi - Intrinsic Impedance
% di - Induction parameter
% ei - Exponential Factor
% ri - Reflection coeficient
% re - Earth R.C.
        
%Step 1 : Calculate basement impedance  
Zn = sqrt(sqrt(-1)*w*mu*resistivities(n)); 
impedances(n) = Zn; 

%Iterate through layers starting from layer j=n-1 (i.e. the layer above the basement)        
for j = n-1:-1:1
    resistivity = resistivities(j);
    thickness = thicknesses(j);
                
    % 3. Compute apparent resistivity from top layer impedance
    %Step 2. Iterate from bottom layer to top(not the basement) 
    % Step 2.1 Calculate the intrinsic impedance of current layer
    dj = sqrt(sqrt(-1)* (w * mu * (1/resistivity)));
    wj = dj * resistivity;
        
    % Step 2.2 Calculate Exponential factor from intrinsic impedance
    ej = exp(-2*thickness*dj);                     

    % Step 2.3 Calculate reflection coeficient using current layer
    %          intrinsic impedance and the below layer impedance
    belowImpedance = impedances(j + 1);
    rj = (wj - belowImpedance)/(wj + belowImpedance); 
    re = rj*ej; 
    Zj = wj * ((1 - re)/(1 + re));
    impedances(j) = Zj;               
end
% Step 3. Compute apparent resistivity from top layer impedance
Z = impedances(1);
absZ = abs(Z); 
apparentResistivity(k) = (absZ * absZ)/(mu * w);
phase(k) = rad2deg(atan2(imag(Z),real(Z)));
end
    