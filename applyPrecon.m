function [ z ] = applyPrecon( Aplusi, Aplusj, rhs, z_in, precon, nx, ny )
%APPLYPRECON Summary of this function goes here
%   Detailed explanation goes here

z = z_in;

idx = 1;
for y = 1:ny
    for x = 1:nx
        t = rhs(idx);
        
        if x > 1
            t = t - Aplusi(idx - 1) * precon(idx - 1) * z(idx - 1);
        end
        if y > 1
            t = t - Aplusj(idx - nx) * precon(idx - nx) * z(idx - nx);
        end
        
        z(idx) = t * precon(idx);
        
        idx = idx + 1;
        
    end
end

for y = ny:-1:1
    for x = nx:-1:1
        idx = getIdx(x,y,nx);
        t = z(idx);
        
        if x < nx
            t = t - Aplusi(idx) * precon(idx) * z(idx + 1);
        end
        if y < ny
            t = t - Aplusj(idx) * precon(idx) * z(idx + nx);
        end
        
        z(idx) = t * precon(idx);
    end
end




end

