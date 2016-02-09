function [ z ] = applyPrecon( Aplusi, Aplusj, rhs, precon, nx, ny )
%APPLYPRECON Summary of this function goes here
%   Detailed explanation goes here

q = zeros(size(precon,1),1);
z = zeros(size(precon,1),1);

idx = 1;
for y = 1:ny
    for x = 1:nx
        t = rhs(idx);
        
        if x > 1
            t = t - Aplusi(idx - 1) * precon(idx - 1) * q(idx - 1);
        end
        if y > 1
            t = t - Aplusj(idx - nx) * precon(idx - nx) * q(idx - nx);
        end
        
        q(idx) = t * precon(idx);
        
        idx = idx + 1;
        
    end
end

for y = ny-1:-1:1
    for x = nx-1:-1:1
        idx = getIdx(x,y,nx);
        t = q(idx);
        
        if x < nx-1
            t = t - Aplusi(idx) * precon(idx) * z(idx + 1);
        end
        if y < ny-1
            t = t - Aplusj(idx) * precon(idx) * z(idx + nx);
        end
        
        z(idx) = t * precon(idx);
    end
end




end

