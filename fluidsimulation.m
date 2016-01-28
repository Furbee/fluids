%% velocity field





%% declaration

g = 9.82; % gravity
rho = 1e3; % density (1e3 for water, 1.3 for air)
dt = 1e-2; % time step
tf = 4; % final time
nx = 90; % number of x-gridpoints
ny = 90; % number of y-gridpoints
lxy = 1; % size of each grid

%% create grid

p = zeros(nx,ny); % pressure at each grid
u = zeros(nx+1,ny); % speed in x-direction
v = zeros(nx,ny+1); % speed in y-direction
pn = zeros(nx,ny); % next pressure at each grid
un = zeros(nx+1,ny); % next speed in x-direction
vn = zeros(nx,ny+1); % next speed in y-direction

%for t=1:10

%% advection

umax = max(max(u),max(v+sqrt(5*lxy*g))); % update max speed
dt = lxy/umax; % update dt
dxy = 0.5;

% calc new positions
for a = 1:nx
    for b = 1:ny
        %pn(a,b) = p(a,b) - dt*(dx*u(a,b)+dx*u(a+1,b));
        %pn(a,b) = pn(a,b) + (p(a,b) - dt*(dx*v(a,b)+dx*u(a,b+1)));
        tx = a - dt*(u(a,b));
        ty = b - dt*(v(a,b));
        alphax = (tx - a);
        alphay = (ty - b);
        pn(a,b) = (
    end
end
p=pn;

%% project

scale = dt/(rho*dx);
for a = 1:nx
    for b = 1:ny
        un(a,b) = u(a,b) - scale * p(a,b);
        un(a,b) = u(a,b) + scale * p(
    end
end

%end
