function [ q11, q21, q12, q22 ] = getValues( tx, ty, w, quantity )
%GETVALUES Summary of this function goes here
%   Detailed explanation goes here



q12 = quantity( getIdx( floor(tx),   floor(ty), w)    );
q22 = quantity( getIdx( ceil(tx),    floor(ty), w)    );
q11 = quantity( getIdx( floor(tx),   ceil(ty), w)     );
q21 = quantity( getIdx( ceil(tx),    ceil(ty), w)     );


end

