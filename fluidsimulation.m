%% declaration

g = -9.82; % gravity
rho = 1.3; % density (1e3 for water, 1.3 for air)
dt = 1e-2; % time step
tf = 4; % final time
nx = 90; % number of x-gridpoints
ny = 90; % number of y-gridpoints
lxy = 1; % size of each grid

%% create grid

p = zeros(nx*ny, 1); % pressure at each grid
d = zeros(nx*ny, 1);
dn = zeros(nx*ny, 1);
u = ones((nx+1)*ny, 1); % speed in x-direction
v = zeros(nx*(ny+1), 1); % speed in y-direction
pn = zeros(nx*ny, 1); % next pressure at each grid
un = zeros((nx+1)*ny, 1); % next speed in x-direction
vn = zeros(nx*(ny+1), 1); % next speed in y-direction
rhs = zeros(nx*ny, 1); % right hand side
Adiag = zeros(nx*ny, 1); % Coeffcient matrix for pressure equations
Aplusi = zeros(nx*ny, 1); %
Aplusj = zeros(nx*ny, 1); %

% MIC
tau_mic = 0.97;
sigma_mic = 0.25;
precon = zeros(nx*ny,1);

% Apply precon
q = zeros(nx*ny,1);
z = zeros(nx*ny,1);

s = zeros(nx*ny,1); % Search vector (matrix)
iter_limit = 10;


colormap winter

for outer_t=1:40
    
    
    umax = max(max(max(u)),max(max(v+sqrt(5*lxy*abs(g))))); % update max speed
    dt = lxy/umax; % update dt
%     dt = 0.0025;
    dxy = 0.5;
    
    %% advect
    
    
    
    % calc new positions
    
    idx = 1;
    for a = 1:nx
        for b = 1:ny
            tx = a - dt*(u(idx));
            ty = b - dt*(v(idx));
            alphax = (tx - floor(tx));
            alphay = (ty - floor(ty));
            un(idx) = (1-alphax)*u(idx) + alphax*u(idx + 1);
            vn(idx) = (1-alphay)*v(idx) + alphay*v(idx + nx);            
            idx = idx + 1;
        end
    end
    
    % update u and v
    u = un;
    v = vn;
    
    
    %% forces
    
    idx = 1;
    for a = 1:nx
        for b = 1:ny
            vn(idx) = v(idx) + dt*g;
            idx = idx + 1;
        end
    end
    
    % update u and v
    v = vn;
    
    
    %% project
    
    rhs = zeros(nx*ny,1); % right hand side
    Adiag = zeros(nx*ny,1); % Coeffcient matrix for pressure equations
    Aplusi = zeros(nx*ny,1); %
    Aplusj = zeros(nx*ny,1); %
    
    % Calculate negative divergence (fig 4.2 in Bridson)
    scale = 1/dxy;
    for a = 2:nx-2
        for b = 2:ny-2
            idx = getIdx(a,b,nx);
            rhs(idx) = -scale * ((u(idx + 1) - u(idx)) ...
                + (v(idx + nx) - v(idx)));
        end
    end
    
    % Modify RHS for solid velocities (fig. 4.3 in Bridson)
    %     for a = 2:nx-1
    %         for b = 2:ny-1
    %             if a == 2
    %                 rhs(a,b) = rhs(a,b) - (scale * u(a,b));
    %             end
    %             if a == nx-1
    %                 rhs(a,b) = rhs(a,b) + (scale * u(a+1,b));
    %             end
    %             if b == 2
    %                 rhs(a,b) = rhs(a,b) - (scale * v(a,b));
    %             end
    %             if b == ny-1
    %                 rhs(a,b) = rhs(a,b) + (scale * v(a,b+1));
    %             end
    %         end
    %     end
    
    % Set up matrix entities for the pressure equations
    scale = dt / (rho * dxy * dxy);
    
    % i
    for a = 2:nx-2
        for b = 2:ny-1
            idx = getIdx(a,b,nx);
            Adiag(idx) = Adiag(idx) + scale;
            Adiag(idx+1) = Adiag(idx+1) + scale;
            Aplusi(idx) = -scale;
        end
    end
    
    % j
    for a = 2:nx-1
        for b = 2:ny-2
            idx = getIdx(a,b,nx);
            Adiag(idx) = Adiag(idx) + scale;
            Adiag(idx + nx) = Adiag(idx + nx) + scale;
            Aplusj(idx) = -scale;
        end
    end
    
    
    % Apply precon
    q = zeros(nx*ny,1);
    z = zeros(nx*ny,1);
    
    % MIC(0) preconditioner
    for a = 2:nx-1
        for b = 2:ny-1
            idx = getIdx(a,b,nx);
            e = Adiag(idx);
            
            px = Aplusi(idx - 1) * precon(idx - 1);
            py = Aplusj(idx - 1) * precon(idx - 1);
            e = e - (px*px + tau_mic*px*py);
            
            px = Aplusi(idx - nx) * precon(idx - nx);
            py = Aplusj(idx - nx ) * precon(idx - nx);
            e = e - (py*py + tau_mic*px*py);
            
            if e < sigma_mic * Adiag(idx)
                e = Adiag(idx);
            end
            
            precon(idx) = 1/sqrt(e);
            
        end
    end
    
    
    % Apply precon
    for a = 2:nx-1
        for b = 2:ny-1
            idx = getIdx(a,b,nx);
            t = rhs(idx);
            t = t - Aplusi(idx - 1) * precon(idx - 1) * q(idx - 1);
            t = t - Aplusj(idx - nx) * precon(idx - nx) * q(idx - nx);
            q(idx) = t * precon(idx);
        end
    end
    for a = nx-1:-1:2
        for b = ny-1:-1:2
            idx = getIdx(a,b,nx);
            t = q(idx);
            t = t - Aplusi(idx) * precon(idx) * z(idx + 1);
            t = t - Aplusj(idx) * precon(idx) * z(idx + nx);
            z(idx) = t * precon(idx);
        end
    end
    % End apply precon
    
    
    s = z;
    
    
    maxError = max(norm(rhs));
    if maxError < 1e-5
         continue;
    end
    
    sigma = dot(rhs,z);
    
    % Iterative solver
    for iter = 1:iter_limit
        
        % Matrix vector product
        for a = 2:nx-2
            for b = 2:ny-2
                idx = getIdx(a,b,nx);
                t = Adiag(idx) * s(idx);
                
                t = t + Aplusi(idx - 1) * s(idx - 1);
                t = t + Aplusj(idx - nx) * s(idx - nx);
                t = t + Aplusi(idx) * s(idx + 1);
                t = t + Aplusj(idx) * s(idx + nx);
                
                z(idx) = t;
                
            end
        end
        
        
        alpha = sigma / dot(z,s);
        
        % Scaled add
        for idx = 1:nx*ny
            p(idx) = p(idx) + s(idx) * alpha;
            rhs(idx) = rhs(idx) - z(idx) * alpha;
        end
        
        maxError = max(norm(rhs));
        if maxError < 1e-5
            %printf('Exiting solver after %d iterations, maximum error is %f\n', iter, maxError);
            %break;
        end
        
        % Apply precon
        for a = 2:nx-1
            for b = 2:ny-1
                idx = getIdx(a,b,nx);
                t = rhs(idx);
                t = t - Aplusi(idx - 1) * precon(idx - 1) * q(idx - 1);
                t = t - Aplusj(idx - nx) * precon(idx - nx) * q(idx - nx);
                q(idx) = t * precon(idx);
            end
        end
        for a = nx-1:-1:2
            for b = ny-1:-1:2
                idx = getIdx(a,b,nx);
                t = q(idx);
                t = t - Aplusi(idx) * precon(idx) * z(idx + 1);
                t = t - Aplusj(idx) * precon(idx) * z(idx + nx);
                z(idx) = t * precon(idx);
            end
        end
        % End apply precon
        
        sigmaNew = dot(z, rhs);
        
        % Scaled add
        for idx = 1:nx*ny
            s(idx) = z(idx) + s(idx) * sigmaNew/sigma;
        end
        
        sigma = sigmaNew;
        
    end
    
    
    
    
    
    
    %     A = delsq(numgrid('S',92));
    %
    %     L = ichol(A,struct('michol','on'));
    %     [p, flag1,rr1,iter1,rv1] = pcg(A,rhs(:),0.01,50,L,L');
    %
    %     temp = rhs(:);
    %
    %     p = reshape(p,[90 90]);
    
    
    % Pressure update
    scale = dt/(rho*dxy);
    for a = 2:nx-2
        for b = 2:ny-2
            idx = getIdx(a,b,nx);
            un(idx) = u(idx) - scale * p(idx);
            un(idx + 1) = u(idx + 1) + scale * p(idx);
            vn(idx) = v(idx) - scale * p(idx);
            vn(idx + nx) = v(idx + nx) + scale *p(idx);
        end
    end
    
    % Boundaries, x
    for a = 1:nx+1
        for b = 1:ny
            idx = getIdx(a,1,nx);
            un(idx) = 0;
            idx = getIdx(a,ny,nx);
            un(idx) = 0;
            idx = getIdx(1,b,nx);
            un(idx) = 0;
            idx = getIdx(nx+1,b,nx);
            un(idx) = 0;
        end
    end
    
    % Boundaries, y
    for a = 1:nx
        for b = 1:ny+1
            idx = getIdx(1,b,nx);
            vn(idx) = 0;
            idx = getIdx(nx,b,nx);
            vn(idx) = 0;
            idx = getIdx(a,1,nx);
            vn(idx) = 0;
            idx = getIdx(a,ny+1,nx);
            vn(idx) = 0;
        end
    end
    
    
    u = un;
    v = vn;
    
    
    
    
    
    
    
    %imshowpair(u',v');
    
    temp_p = reshape(p, [nx ny]);
    temp_u = reshape(u, [ny nx+1]);
    
    imagesc(temp_p)
    drawnow
end