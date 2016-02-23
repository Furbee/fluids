clear

%% declaration

g = 9.82; % gravity
rho = 0.1; % density (1e3 for water, 1.3 for air)
densitySoot = 0.11;
densityAir = rho;
tAmb = 273;
kDiss = 0.1;
epsilon = 1;

dt = 1e-2; % time step
tf = 4; % final time
nx = 64; % number of x-gridpoints
ny = 64; % number of y-gridpoints
%lxy = 1; % size of each grid
lxy = 1.0/min(nx,ny);
maxtime = 10.0; % set simulation length
elapsedtime = 0.0;
time = 0.0;
pasttime = 0.0;

%% create grid

p = zeros(nx*ny, 1); % pressure at each grid
d = zeros(nx*ny, 1); % soot concentration
dn = zeros(nx*ny, 1);
T = zeros(nx*ny, 1); % temperature
Tn = zeros(nx*ny, 1);
u = zeros((nx+1)*ny, 1); % speed in x-direction
v = zeros(nx*(ny+1), 1); % speed in y-direction
pn = zeros(nx*ny, 1); % next pressure at each grid
un = zeros((nx+1)*ny, 1); % next speed in x-direction
vn = zeros(nx*(ny+1), 1); % next speed in y-direction
w = zeros(nx*ny, 1); % vorticity
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



disp('Writing to newfile.avi...')
video = VideoWriter('newfile.avi','Uncompressed AVI');
open(video)

while time < maxtime
    tic;
    %disp('NEW LOOP')
    %umax = max(max(max(u)),max(max(v+sqrt(5*lxy*abs(g))))); % update max speed
    %dt = lxy/umax; % update dt
    dt = 1/60;
    %dxy = 0.5;
    
    
    
    [ d ] = addInFlow( 0.45, 0.80, 0.60, 0.83, nx, ny, 0.5, 0.5, ...
        dxy, 1.0, d);
    
    [ T ] = addInFlow( 0.45, 0.80, 0.60, 0.83, nx, ny, 0.5, 0.5, ...
        dxy, tAmb + 300.0, T);
    
    [ u ] = addInFlow( 0.45, 0.80, 0.60, 0.83, nx+1, ny, 0.0, 0.5, ...
        dxy, 0.0, u);
    [ v ] = addInFlow( 0.45, 0.80, 0.60, 0.83, nx, ny+1, 0.5, 0.0, ...
        dxy, 0.0, v);
    
    
    
    %% external forces
    
%         for y = 1:ny+1
%             for x = 1:nx
%                 idx = getIdx(x,y,nx);
%                 v(idx) = v(idx) + dt*g;
%             end
%         end

    % buoyancy
    
    alpha = (densitySoot - densityAir)/densityAir;
    
    idx = 1;
    for y = 1:ny
        for x = 1:nx
            
            buoyancy = dt * g * (alpha*d(idx)- (T(idx) - tAmb)/tAmb);
            
            v(idx) = v(idx) + buoyancy*0.5;
            v(idx+nx) = v(idx+nx) + buoyancy*0.5;
            
            idx = idx + 1;
            
        end
    end
    
    
    for x = 1:nx
        idx = getIdx(x,1,nx);
        v(idx) = 0.0;
        v(idx+1) = 0.0;
        idx = getIdx(x,ny,nx);
        v(idx) = 0.0;
        v(idx+nx) = 0.0;
    end
    
    
    % vorticity confinement
    
    % for every cell
    % 1. get velocities in surrounding cell centers
    % 2. calculate vorticity, w, in each cell using central differences
    % 3. calculate gradient of the vorticity in each cell using central
    %    differences
    % 4. calculate unit vector N by normalizing w
    % 5. calculate force with F = epsilon * dxy * ( N x w)
    % 6. apply force to every cell face
    
%    for y = 2:ny-1
%        for x = 2:nx-1
%            uRight = (u(getIdx(x+2,y,nx+1)) + u(getIdx(x+1,y,nx+1)))*0.5;
%            uLeft = (u(getIdx(x,y,nx+1)) + u(getIdx(x-1,y,nx+1)))*0.5;
%            vTop = (v(getIdx(x,y-1,nx)) + v(getIdx(x,y,nx)))*0.5;
%            vBottom = (v(getIdx(x,y+2,nx)) + v(getIdx(x,y+1,nx)))*0.5;
%            
%            w(getIdx(x,y,nx)) = ( (uRight - uLeft) - (vTop - vBottom)) / ( 2 * dxy);
%            
%        end
%    end
%    
%    for y = 3:ny-2
%        for x = 3:nx-2
%            
%            wi = w(getIdx(x+1,y,nx)) - w(getIdx(x-1,y,nx));
%            wj = w(getIdx(x,y+1,nx)) - w(getIdx(x,y-1,nx));
%            gradW = [ wi wj ];
%            gradW = gradW / (2 * dxy);
%            M = 1 / (dxy * dt);
%            
%            N = gradW / ( norm(gradW) + 10e-20 * M);
%            
%            fconf = epsilon * dxy * cross(N, gradW);
%            
%            
%            
%            
%            
%        end
%    end
%    
%    
    
    %% project
    
    
    % Calculate negative divergence (fig 4.2 in Bridson)
    scale = 1/dxy;
    
    idx = 1;
    for y = 1:ny
        for x = 1:nx
            %             idx = getIdx(x,y,nx);
            rhs(idx) = -scale * ((u(getIdx(x+1,y,nx+1)) - u(getIdx(x,y,nx+1))) ...
                + (v(getIdx(x,y+1,ny)) - v(getIdx(x,y,ny))));
            assert(isnan(rhs(idx)) == 0)
            idx = idx + 1;
        end
    end
    %rhs_max = max(abs(rhs))
    
    
       %  [ p ] = project2( rhs, nx, ny, dt, rho, dxy, iter_limit);
    
    
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
    Adiag = zeros(nx*ny, 1); % Coeffcient matrix for pressure equations
    %     Aplusi = zeros(nx*ny, 1); %
    %     Aplusj = zeros(nx*ny, 1); %
    
    idx = 1;
    for y = 1:ny
        for x = 1:nx
            if x < nx
                Adiag(idx) = Adiag(idx) + scale;
                Adiag(idx+1) = Adiag(idx+1) + scale;
                Aplusi(idx) = (-scale);
            else
                Aplusi(idx) = 0.0;
            end
            
            if y < ny
                Adiag(idx) = Adiag(idx) + scale;
                Adiag(idx + nx) = Adiag(idx + nx) + scale;
                Aplusj(idx) = (-scale);
            else
                Aplusj(idx) = 0.0;
            end
            
            %idx = getIdx(x,y,nx);
            idx = idx + 1;
            
        end
    end
    %Adiag_max = max(abs(Adiag))
    %Aplusi_max = max(abs(Aplusi))
    %Aplusj_max = max(abs(Aplusj))
    
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
    
    %precon_max = max(abs(precon))
    
    [p, rhs] = project( Adiag, Aplusi, Aplusj, rhs,  precon, nx, ny, iter_limit );
    %p_max = max(abs(p))

    % Pressure update
    scale = dt/(rho*dxy);
    idx = 1;
    for y = 1:ny
        for x = 1:nx
            uvidx = getIdx(x,y,nx+1);
            u(uvidx) = u(uvidx) - scale * p(idx);
            u(uvidx+1) = u(uvidx+1) + scale * p(idx);
            uvidx = getIdx(x,y,nx);
            v(uvidx) = v(uvidx) - scale * p(idx);
            v(uvidx + nx) = v(uvidx + nx) + scale * p(idx);
            idx = idx + 1;
        end
    end

    %u_max = max(abs(u))
    %v_max = max(abs(v))
    % Boundaries, x
    for y = 1:ny
        idx = getIdx(1,y,nx+1);
        u(idx) =  0.0;
        u(idx+1) = 0.0;
        idx = getIdx(nx,y,nx+1);
        u(idx) = 0.0;
        u(idx+1) = 0.0;
    end
    
    for x = 1:nx
        idx = getIdx(x,1,nx);
        v(idx) = 0.0;
        v(idx+1) = 0.0;
        idx = getIdx(x,ny,nx);
        v(idx) = 0.0;
        v(idx+nx) = 0.0;
    end
    
    
    %% advect
    % calc new positions
    
    
    
    idx = 1;
    for y = 1:ny
        for x = 1:nx
            
            % offset
            ix = x + 0.5;
            iy = y + 0.5;
            
            [x0, y0] = rungeKutta3( ix, iy, dt, u, v, dxy, nx, ny);
            %dn(idx) = lerp2(x0, y0, 0.5, 0.5, nx, ny, d);
            dn(idx) = cerp2(x0, y0, nx, ny, 0.5, 0.5, d);
            
            dn(idx) = max(0, dn(idx) * exp( -kDiss * dt));
            
            Tn(idx) = cerp2(x0, y0, nx, ny, 0.5, 0.5, T);
            
            idx = idx + 1;
        end
    end
    
    
    idx = 1;
    for y = 1:ny
        for x = 1:nx+1
            
            % offset
            ix = x + 0.0;
            iy = y + 0.5;
            
            [x0, y0] = rungeKutta3( ix, iy, dt, u, v, dxy, nx, ny);

            %un(idx) = lerp2(x0, y0, 0.0, 0.5, nx+1, ny, u);
            un(idx) = cerp2(x0, y0, nx+1, ny, 0.0, 0.5, u);
            idx = idx + 1;
        end
    end
    
    idx = 1;
    for y = 1:ny+1
        for x = 1:nx
            
            % offset
            ix = x + 0.5;
            iy = y + 0.0;
            
            [x0, y0] = rungeKutta3( ix, iy, dt, u, v, dxy, nx, ny);
            %vn(idx) = lerp2(x0, y0, 0.5, 0.0, nx, ny+1, v);
            vn(idx) = cerp2(x0, y0, nx, ny+1, 0.5, 0.0, v );
            idx = idx+1;
        end
    end
    
    
    % update u and v
    d = dn;
    T = Tn;
    u = un;
    v = vn;
    
    %u_max = max(abs(u))
    %v_max = max(abs(v))
    
    %imshowpair(u',v');
    
    
    
    
    %     u(getIdx(10,10,nx)) = 100;
    %     u(getIdx(10,12,nx)) = 100;
    %
    %     pause(2)
    
    
    
    temp_d = reshape(d, [ny, nx]);
    temp_u = reshape(u, [ny+1, nx]);
    temp_v = reshape(v, [ny, nx+1]);
    
    %imagesc(temp_u)
    colormap gray
    imagesc(temp_d');
    
    drawnow
    F = getframe;
    writeVideo(video,F);
    
%     imagesc(reshape(p, [ny, nx]));
    
    
    time = time + dt;
    
    perc = time/maxtime;
    elapsedtime = elapsedtime + (toc + pasttime);
    tic;
    tottime = (1/perc)*elapsedtime;
    mins = (tottime - elapsedtime)/60;
    fmins = floor(mins);
    secs = floor((mins-fmins)*60);
    str = sprintf('Estimated remaining time: %i minutes and %i seconds', fmins,secs);
    disp(str)
    pasttime = toc;
    
end
close(video)
disp('Done!')





