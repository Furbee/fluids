%% declaration

g = -9.82; % gravity
rho = 1e3; % density (1e3 for water, 1.3 for air)
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
un = ones((nx+1)*ny, 1); % next speed in x-direction
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
iter_limit = 80;


colormap winter

for outer_t=1:100
    
    
    umax = max(max(max(u)),max(max(v+sqrt(5*lxy*abs(g))))); % update max speed
    dt = lxy/umax; % update dt
    %     dt = 0.0025;
    dxy = 0.5;
    
    %% advect
    
    
    
    % calc new positions
    
    idx = 1;
    for y = 1:ny
        for x = 1:nx
            tx = x - dt*(u(idx));
            ty = y - dt*(v(idx));
            
            tx = min(max(tx, 1), nx);
            ty = min(max(ty, 1), ny);
            
            [ q11, q21, q12, q22 ] = getValues( tx, ty, nx, u );
            un(idx) = lerp(tx, ty, q11, q21, q12, q22);
            
            [ q11, q21, q12, q22 ] = getValues( tx, ty, nx, v );
            vn(idx) = lerp(tx, ty, q11, q21, q12, q22 );
            idx = idx + 1;
        end
    end
    
    
    
    % update u and v
    u = un;
    v = vn;
    
    
    %% forces
    
    idx = 1;
    for y = 1:ny+1
        for x = 1:nx
            vn(idx) = v(idx) + dt*g;
            idx = idx + 1;
        end
    end
    
    % update u and v
    v = vn;
    
    
    %% project
    
    
    % Calculate negative divergence (fig 4.2 in Bridson)
    scale = 1/dxy;
    for y = 1:ny
        for x = 1:nx
            idx = getIdx(x,y,nx);
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
    
    idx = 1;
    for y = 1:ny
        for x = 1:nx
            if x < nx
                Adiag(idx) = Adiag(idx) + scale;
                Adiag(idx+1) = Adiag(idx+1) + scale;
                Aplusi(idx) = (-scale);
            end
            
            if y < ny
                Adiag(idx) = Adiag(idx) + scale;
                Adiag(idx + nx) = Adiag(idx + nx) + scale;
                Aplusj(idx) = (-scale);
            end
            
            %idx = getIdx(x,y,nx);
            idx = idx + 1;
            
        end
    end
    
    % MIC(0) preconditioner
    idx = 1;
    for y = 1:ny
        for x = 1:nx
            %             idx = getIdx(x,y,nx);
            e = Adiag(idx);
            
            if x > 1
                px = Aplusi(idx - 1) * precon(idx - 1);
                py = Aplusj(idx - 1) * precon(idx - 1);
                e = e - (px*px + tau_mic*px*py);
            end
            
            if y > 1
                px = Aplusi(idx - nx) * precon(idx - nx);
                py = Aplusj(idx - nx ) * precon(idx - nx);
                e = e - (py*py + tau_mic*px*py);
            end
            
            if e < sigma_mic * Adiag(idx)
                e = Adiag(idx);
            end
            
            precon(idx) = 1/sqrt(e);
            
            idx = idx + 1;
            
        end
    end
    
    p = project( Adiag, Aplusi, Aplusj, rhs, precon, nx, ny, iter_limit );
    
    % Pressure update
    scale = dt/(rho*dxy);
    for y = 1:ny
        for x = 1:nx
            idx = getIdx(x,y,nx);
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
    
    temp_u = reshape(u, [ny nx+1]);
    
    imagesc(temp_u)
    drawnow
end