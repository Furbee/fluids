function [ xy ] = lerp( a, b, x )
%LERP Summary of this function goes here
%   Detailed explanation goes here

xy = a*(1.0 - x) + b*x;


end

