function [ x_out, y_out ] = rungeKutta3( x, y, dt, u, v, dxy, w, h )
%RK3 Summary of this function goes here
%   Detailed explanation goes here

firstU = lerp2(x, y, 0.0, 0.5, w+1, h, u)/dxy;
firstV = lerp2(x, y, 0.5, 0.0, w, h+1 ,v)/dxy;

midX = x - 0.5*dt*firstU;
midY = y - 0.5*dt*firstV;

midU = lerp2(midX, midY, 0.0, 0.5, w+1, h, u)/dxy;
midV = lerp2(midX, midY, 0.5, 0.0, w, h+1, v)/dxy;

lastX = x - 0.75*dt*midU;
lastY = y - 0.75*dt*midV;

lastU = lerp2(lastX, lastY, 0.0, 0.5, w+1, h, u);
lastV = lerp2(lastX, lastY, 0.5, 0.0, w, h+1, v);

x_out = x - dt*((2.0/9.0)*firstU + (3.0/9.0)*midU + (4.0/9.0)*lastU);
y_out = y - dt*((2.0/9.0)*firstV + (3.0/9.0)*midV + (4.0/9.0)*lastV);


end

