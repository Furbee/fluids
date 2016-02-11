%% declaration

g = -9.82; % gravity
rho = 0.1; % density (1e3 for water, 1.3 for air)
dt = 1e-2; % time step
tf = 4; % final time
nx = 64; % number of x-gridpoints
ny = 64; % number of y-gridpoints
%lxy = 1; % size of each grid
lxy = 1.0/min(nx,ny);

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
iter_limit = 600;


%colormap winter

%dxy = lxy;
dxy = lxy;



for outer_t=1:100
    
    
    %umax = max(max(max(u)),max(max(v+sqrt(5*lxy*abs(g))))); % update max speed
    %dt = lxy/umax; % update dt
    dt = 0.005;
    %dxy = 0.5;
    
    [ d ] = addInFlow( 0.45, 0.2, 0.55, 0.21, nx, ny, 0.5, 0.5, ...
        dxy, 1.0, d);
    
    [ u ] = addInFlow( 0.45, 0.2, 0.55, 0.21, nx, ny, 0.0, 0.5, ...
        dxy, 0.0, u);
    [ v ] = addInFlow( 0.45, 0.2, 0.55, 0.21, nx, ny, 0.5, 0.0, ...
        dxy, 3.0, v);
    
    
    
    %% forces
    
    %     for y = ny/4:ny/2
    %         for x = nx/4:nx/2
    %             idx = getIdx(x,y,nx);
    %             v(idx) = v(idx) + dt*g;
    %         end
    %     end
    %
    %
    %     % update u and v
    %     v = vn;
    
    
    %% project
    
    
    % Calculate negative divergence (fig 4.2 in Bridson)
    scale = 1.0/dxy;
    
    idx = 1;
    for y = 1:ny
        for x = 1:nx
            %             idx = getIdx(x,y,nx);
            rhs(idx) = -scale * ((u(getIdx(x+1,y,nx)) - u(getIdx(x,y,nx))) ...
                + (v(getIdx(x,y+1,nx)) - v(getIdx(x,y,nx))));
            assert(isnan(rhs(idx)) == 0)
            
            idx = idx + 1;
        end
    end
    
    
         [ p ] = project2( rhs, nx, ny, dt, rho, dxy, iter_limit);
    
%     
%     % Modify RHS for solid velocities (fig. 4.3 in Bridson)
%     %     for a = 2:nx-1
%     %         for b = 2:ny-1
%     %             if a == 2
%     %                 rhs(a,b) = rhs(a,b) - (scale * u(a,b));
%     %             end
%     %             if a == nx-1
%     %                 rhs(a,b) = rhs(a,b) + (scale * u(a+1,b));
%     %             end
%     %             if b == 2
%     %                 rhs(a,b) = rhs(a,b) - (scale * v(a,b));
%     %             end
%     %             if b == ny-1
%     %                 rhs(a,b) = rhs(a,b) + (scale * v(a,b+1));
%     %             end
%     %         end
%     %     end
%     
%     % Set up matrix entities for the pressure equations
%     scale = dt / (rho * dxy * dxy);
%     Adiag = zeros(nx*ny, 1); % Coeffcient matrix for pressure equations
%     %     Aplusi = zeros(nx*ny, 1); %
%     %     Aplusj = zeros(nx*ny, 1); %
%     
%     idx = 1;
%     for y = 1:ny
%         for x = 1:nx
%             if x < nx
%                 Adiag(idx) = Adiag(idx) + scale;
%                 Adiag(idx+1) = Adiag(idx+1) + scale;
%                 Aplusi(idx) = (-scale);
%             else
%                 Aplusi(idx) = 0.0;
%             end
%             
%             if y < ny
%                 Adiag(idx) = Adiag(idx) + scale;
%                 Adiag(idx + nx) = Adiag(idx + nx) + scale;
%                 Aplusj(idx) = (-scale);
%             else
%                 Aplusj(idx) = 0.0;
%             end
%             
%             %idx = getIdx(x,y,nx);
%             idx = idx + 1;
%             
%         end
%     end
%     
%     % MIC(0) preconditioner
%     idx = 1;
%     for y = 1:ny
%         for x = 1:nx
%             %             idx = getIdx(x,y,nx);
%             e = Adiag(idx);
%             
%             if x > 1
%                 px = Aplusi(idx - 1) * precon(idx - 1);
%                 py = Aplusj(idx - 1) * precon(idx - 1);
%                 e = e - (px*px + tau_mic*px*py);
%             end
%             
%             if y > 1
%                 px = Aplusi(idx - nx) * precon(idx - nx);
%                 py = Aplusj(idx - nx ) * precon(idx - nx);
%                 e = e - (py*py + tau_mic*px*py);
%             end
%             
%             if e < sigma_mic * Adiag(idx)
%                 e = Adiag(idx);
%             end
%             
%             precon(idx) = 1/sqrt(e);
%             
%             idx = idx + 1;
%             
%         end
%     end
%     
%     [p, rhs] = project( Adiag, Aplusi, Aplusj, rhs,  precon, nx, ny, iter_limit );
    
    % Pressure update
    scale = dt/(rho*dxy);
    idx = 1;
    for y = 1:ny
        for x = 1:nx
            %idx = getIdx(x,y,nx);
            uvidx = getIdx(x,y,nx);
            u(uvidx) = u(uvidx) - scale * p(idx);
            u(uvidx+1) = u(uvidx+1) + scale * p(idx);
            v(uvidx) = v(uvidx) - scale * p(idx);
            v(uvidx + nx) = v(uvidx + nx) + scale * p(idx);
            idx = idx + 1;
        end
    end
    
    % Boundaries, x
    for y = 1:ny
        idx = getIdx(1,y,nx);
        u(idx) =  0.0;
        idx = getIdx(nx,y,nx);
        u(idx) = 0.0;
    end
    
    for x = 1:nx
        idx = getIdx(x,1,nx);
        v(idx) = 0.0;
        idx = getIdx(x,ny,nx);
        v(idx) = 0.0;
    end
    
    
    %% advect
    % calc new positions
    
    idx = 1;
    for y = 1:ny
        for x = 1:nx
            
            % offset
            ix = x + 0.5;
            iy = y + 0.5;
            
            [x0, y0] = rungeKutta3( ix, iy, 0.5, 0.5, dt, u, v, dxy, nx, ny);
            dn(idx) = cerp2(x0, y0, nx, ny, 0.5, 0.5, d);
            idx = idx + 1;
        end
    end
    
    idx = 1;
    for y = 1:ny
        for x = 1:nx+1
            
            % offset
            ix = x + 0.0;
            iy = y + 0.5;
            
            [x0, y0] = rungeKutta3( ix, iy, 0.0, 0.5, dt, u, v, dxy, nx+1, ny);
            un(idx) = cerp2(x0, y0, nx+1, ny, 0.0, 0.5, u);
            idx = idx + 1;
        end
    end
    
    idx = 1;
    for y = 1:ny+1
        for x = 1:nx
            
            % offset
            ix = ix + 0.5;
            iy = iy + 0.0;
            
            [x0, y0] = rungeKutta3( ix, iy, 0.5, 0.0, dt, u, v, dxy, nx, ny+1);
            vn(idx) = cerp2(x0, y0, nx, ny+1, 0.5, 0.0, v );
            idx = idx+1;
        end
    end
    
    % update u and v
    d = dn;
    u = un;
    v = vn;
    
    %imshowpair(u',v');
    
    temp_d = reshape(d, [ny, nx]);
    temp_u = reshape(u, [ny, nx+1]);
    temp_v = reshape(v, [nx, ny+1]);
    
    %     imagesc(temp_u)
    imagesc(temp_v');
    
    %imagesc(reshape(p, [ny, nx]));
    
    drawnow
    
    
end







