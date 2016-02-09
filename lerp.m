function [ xy ] = lerp( tx, ty, q11, q21, q12, q22 )
%LERP Summary of this function goes here
%   Detailed explanation goes here


alphax = (tx - floor(tx));
alphay = (ty - floor(ty));

xy1 = (1 - alphax)  * q11 + alphax * q21;
xy2 = (1 - alphax)  * q12 + alphax * q22;
xy  = (1 - alphay)  * xy1 + alphay * xy2;



end

