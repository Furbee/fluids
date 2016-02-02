function [ idx ] = getIdx( i, j, width )
%GETIDX Get index of specified cell

    idx = i + (j-1)*width;

end


% 
%   1   4   7
%   2   5   8
%   3   6   9
%
%
% getIdx(2,2,3) = 2 + (2-1)*3 = 5 (yep)
% getIdx(2,3,3) = 2 + (3-1)*3 = 8
% getIdx(3,3,3) = 3 + (3-1)*3 = 3 + 6 = 9