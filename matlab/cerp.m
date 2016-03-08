function [ xy ] = cerp( a, b, c, d, x )
%CERP Summary of this function goes here
%   Detailed explanation goes here

xsq = x*x;
xcu = xsq*x;
minV = min(a, min(b, min(c, d)));
maxV = max(a, max(b, max(c, d)));

t = ...
a*(0.0 - 0.5*x + 1.0*xsq - 0.5*xcu) + ...
b*(1.0 + 0.0*x - 2.5*xsq + 1.5*xcu) + ...
c*(0.0 + 0.5*x + 2.0*xsq - 1.5*xcu) + ...
d*(0.0 + 0.0*x - 0.5*xsq + 0.5*xcu);

xy = min(max(t, minV), maxV);


end

