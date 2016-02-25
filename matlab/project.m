function [ p, rhs ] = project( Adiag, Aplusi, Aplusj, rhs, precon, nx, ny, ...
    iter_limit)
%PROJECT Summary of this function goes here
%   Detailed explanation goes here

p = zeros(size(rhs));
z = zeros(size(rhs));

% Apply precon

z = applyPrecon(Aplusi, Aplusj, rhs, z, precon, nx, ny);

s = z;


maxError = max(norm(rhs));
if maxError < 1e-5
    return;
end

sigma = dot(rhs,z)

% Iterative solver
for iter = 1:iter_limit
    
    %         dbstop if any(isnan(Adiag)) == 1
    %         dbstop if any(isnan(Aplusi)) == 1
    %         dbstop if any(isnan(Aplusj)) == 1
    %         dbstop if any(isnan(s)) == 1
    
    % Matrix vector product
    idx = 1;
    for y = 1:ny
        for x = 1:nx
            %idx = getIdx(x,y,nx);
            t = Adiag(idx) * s(idx);
            
            if x > 1
                t = t + Aplusi(idx - 1) * s(idx - 1);
            end
            if y > 1
                t = t + Aplusj(idx - nx) * s(idx - nx);
            end
            if x < nx
                t = t + Aplusi(idx) * s(idx + 1);
            end
            if y < ny
                t = t + Aplusj(idx) * s(idx + nx);
            end
            
            z(idx) = t;
            
            idx = idx + 1;
        end
    end
    
    alpha = sigma / dot(z,s);
    
    % Scaled add
    p = p + s*alpha;
    rhs = rhs - z*alpha;
    
    maxError = max(norm(rhs));
    if maxError < 1e-5
        fprintf('Exiting solver after %d iterations, maximum error is %f\n', iter, maxError);
        %disp('Exiting solver after ' << iter + ' iterations');
        return;
    end
    
    % Apply precon
    z = applyPrecon(Aplusi, Aplusj, rhs, z, precon, nx, ny);
    
    sigmaNew = dot(rhs,z);
    
    % Scaled add
    %         for idx = 1:nx*ny
    %             s(idx) = z(idx) + s(idx) * sigmaNew/sigma;
    %         end
    s = z + s * (sigmaNew/sigma);
    
    sigma = sigmaNew;
    
    
    
end

fprintf('Exceeded budget of %d iterations, maximum error was %f\n', iter_limit, maxError);




end
