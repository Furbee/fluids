function [ src_out ] = addInFlow( x0, y0, x1, y1, w, h, ox, oy, ...
    dxy, v, src)
%ADDINFLOW Summary of this function goes here
%   Detailed explanation goes here

src_out = src;


ix0 = floor(x0/dxy - ox);
iy0 = floor(y0/dxy - oy);
ix1 = floor(x1/dxy - ox);
iy1 = floor(y1/dxy - oy);


for y = iy0:iy1
    for x = ix0:ix1
        idx = getIdx(x,y,w);
        if abs(src_out(idx)) < abs(v)
            src_out(idx) = v;
        end
    end
end


end

