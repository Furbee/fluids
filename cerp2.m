function [ xy ] = cerp2(x, y, w, h, ox, oy, quantity )
%CERP2 Summary of this function goes here
%   Detailed explanation goes here

x = min(max(x - ox, 1.0), w - 0.001);
y = min(max(y - oy, 1.0), h - 0.001);
ix = floor(x);
iy = floor(y);
x = x - ix;
y = y - iy;

x0 = max(ix - 1, 1);
x1 = ix;
x2 = ix + 1;
x3 = min(ix + 2, w - 1);

y0 = max(iy - 1, 1);
y1 = iy;
y2 = iy + 1;
y3 = min(iy + 2, h - 1);

q0 = cerp(quantity(getIdx(x0, y0, w)), quantity(getIdx(x1, y0, w)), ...
    quantity(getIdx(x2, y0, w)), quantity(getIdx(x3, y0, w)), x);
q1 = cerp(quantity(getIdx(x0, y1, w)), quantity(getIdx(x1, y1, w)), ...
    quantity(getIdx(x2, y1, w)), quantity(getIdx(x3, y1, w)), x);
q2 = cerp(quantity(getIdx(x0, y2, w)), quantity(getIdx(x1, y2, w)), ...
    quantity(getIdx(x2, y2, w)), quantity(getIdx(x3, y2, w)), x);
q3 = cerp(quantity(getIdx(x0, y3, w)), quantity(getIdx(x1, y3, w)), ...
    quantity(getIdx(x2, y3, w)), quantity(getIdx(x3, y3, w)), x);

xy = cerp(q0, q1, q2, q3, y);


end

