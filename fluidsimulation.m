%% declaration

g = -9.82; % gravity
rho = 1e3; % density (1e3 for water, 1.3 for air)
dt = 1e-2; % time step
tf = 4; % final time
nx = 90; % number of x-gridpoints
ny = 90; % number of y-gridpoints
lxy = 1; % size of each grid

%% create grid

p = zeros(nx,ny); % pressure at each grid
u = ones(nx+1,ny); % speed in x-direction
v = zeros(nx,ny+1); % speed in y-direction
pn = zeros(nx,ny); % next pressure at each grid
un = zeros(nx+1,ny); % next speed in x-direction
vn = zeros(nx,ny+1); % next speed in y-direction
rhs = zeros(nx,ny); % right hand side
Adiag = zeros(nx,ny); % Coeffcient matrix for pressure equations
Aplusi = zeros(nx,ny); %
Aplusj = zeros(nx,ny); %
tau_mic = 0.97;
sigma_mic = 0.25;


for t=1:10
    
    %% advect
    
    umax = max(max(max(u)),max(max(v+sqrt(5*lxy*abs(g))))); % update max speed
    dt = lxy/umax; % update dt
    dxy = 0.5;
    
    % calc new positions
    for a = 1:nx
        for b = 1:ny
            tx = a - dt*(u(a,b));
            ty = b - dt*(v(a,b));
            alphax = (tx - floor(tx));
            alphay = (ty - floor(ty));
            un(a,b) = (1-alphax)*u(a,b) + alphax*u(a+1,b);
            vn(a,b) = (1-alphay)*v(a,b) + alphay*v(a,b+1);
        end
    end
    
    % update u and v
    u = un;
    v = vn;
    
    
    %% forces
    
    for a = 1:nx
        for b = 1:ny
            vn(a,b) = v(a,b) + dt*g;
        end
    end
    
    % update u and v
    v = vn;
    
    
    %% project
    
    % Calculate negative divergence (fig 4.2 in Bridson)
    scale = 1/dxy;
    for a = 2:nx-2
        for b = 2:ny-2
            rhs(a,b) = -scale * ((u(a+1,b) - u(a,b)) ...
                + (v(a,b+1) - v(a,b)));
        end
    end
    
    
    % Modify RHS for solid velocities (fig. 4.3 in Bridson)
    for a = 2:nx-1
        for b = 2:ny-1
            if a == 2
                rhs(a,b) = rhs(a,b) - (scale * u(a,b));
            end
            if a == nx-1
                rhs(a,b) = rhs(a,b) + (scale * u(a+1,b));
            end
            if b == 2
                rhs(a,b) = rhs(a,b) - (scale * v(a,b));
            end
            if b == ny-1
                rhs(a,b) = rhs(a,b) + (scale * v(a,b+1));
            end
        end
    end
    
    % Set up matrix entities for the pressure equations
    scale = dt / (rho * dxy * dxy);
    
    % i
    for a = 2:nx-2
        for b = 2:ny-1
            Adiag(a,b) = Adiag(a,b) + scale;
            Adiag(a+1,b) = Adiag(a+1,b) + scale;
            Aplusi(a,b) = -scale;
        end
    end
    
    % j
    for a = 2:nx-1
        for b = 2:ny-2
            Adiag(a,b) = Adiag(a,b) + scale;
            Adiag(a,b+1) = Adiag(a,b+1) + scale;
            Aplusj(a,b) = -scale;
        end
    end
    
   
    % MIC(0) preconditioner
    
    A = sparse(Adiag);
    
    L = ichol(A,struct('michol','on'));
    [p, flag1,rr1,iter1,rv1] = pcg(A,rhs,0.01,50,L,L');
    
    
    % Pressure update
    scale = dt/(rho*dxy);
    for a = 2:nx-2
        for b = 2:ny-2
            un(a,b) = u(a,b) - scale * p(a,b);
            un(a+1,b) = u(a+1,b) + scale * p(a,b);
            vn(a,b) = v(a,b) - scale * p(a,b);
            vn(a,b+1) = v(a,b+1) + scale *p(a,b);
        end
    end
    
    % Boundaries, x
    for a = 1:nx+1
        for b = 1:ny
            un(a,1) = 0;
            un(a,ny) = 0;
            un(1,b) = 0;
            un(nx+1,b) = 0;
        end
    end
    
    % Boundaries, y
    for a = 1:nx
        for b = 1:ny+1
            vn(1,b) = 0;
            vn(nx,b) = 0;
            vn(a,1) = 0;
            vn(a,ny+1) = 0;
        end
    end
    
    
    
    u = un;
    v = vn;
    
    
    
    
    
    
    %imagesc(v');
    %drawnow
end