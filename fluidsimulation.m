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
    
    scale = dt/(rho*dxy);
    for a = 2:nx-2
        for b = 2:ny-2
            un(a,b) = u(a,b) - scale * p(a,b);
            un(a+1,b) = u(a+1,b) + scale * p(a,b);
            vn(a,b) = v(a,b) - scale * p(a,b);
            vn(a,b+1) = v(a,b+1) + scale *p(a,b);
        end
    end
    for a = 1:nx+1
        for b = 1:ny
            un(a,1) = 0;
            un(a,ny) = 0;
            un(1,b) = 0;
            un(nx+1,b) = 0;
        end
    end
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
    
    scale = 1/dxy;
    for a = 2:nx-2
        for b = 2:ny-2
            rhs(a,b) = -scale * ((u(a+1,b) - u(a,b)) ...
                + (v(a,b+1) - v(a,b)));
        end
    end
    
    for a = 2:nx
        for b = 2:ny
            if a==2 || a == nx
                rhs(a,b) = rhs(a,b) - (scale * u(a,b));
                rhs(a,b) = rhs(a,b) + (scale * u(a+1,b));
            end if b==2 || b == ny
                rhs(a,b) = rhs(a,b) - (scale * v(a,d));
                rhs(a,b) = rhs(a,b) + (scale * v(a,b+1));
            end
        end
    end
    for a = 1:nx+1
        for b = 1:ny
            un(a,1) = 0;
            un(a,ny) = 0;
            un(1,b) = 0;
            un(nx+1,b) = 0;
        end
    end
    
    %imagesc(v');
    %drawnow
end