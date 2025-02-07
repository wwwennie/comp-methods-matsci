% MD code for 2D system of surface atoms using Verlet Algorithm
%  input:  n = number of atoms
%          ein = initial energy
%          acon = parameter in surface potential; strength of interaction
%          nsteps = number of time steps
%          dt = time step
%  output:  un = potential energy at each time step
%           kn = kinetic energy at each time step
%           en = total energy at each time step
%           tn = temperature at each time step
%           pn = pressure at each time step

% NOTE: technically only works for n = 1 atoms at the moment

function[phin,kn,en,xn,yn]= MDSurf(n,ein,acon,nsteps,dt)

% initialize positions and velocities
[x, y, vx, vy] = initsurf(ein, acon);

% calculate initial energy and forces
[phi]= phisurf(acon,x,y);
[fx,fy]= fsurf(acon,x,y);

% now start the time stepping with the verlet algorithm
% initialize variables
xold = zeros(n,1);
yold = zeros(n,1);

xnew = zeros(n,1);
ynew = zeros(n,1);

% first find the positions at t-dt
if n > 1
    for i=1:n
        xold(i) = x(i) - vx(i)*dt + .5*fx(i)*dt^2;
        yold(i) = y(i) - vy(i)*dt + .5*fy(i)*dt^2;
    end
else
    xold = x - vx*dt + .5*fx*dt^2;
    yold = y - vy*dt + .5*fy*dt^2;
end

% start the time steps
for j=1:nsteps
    k = 0;
 %  find positions for time t + dt
 %  find velocities for time t
 %  find kinetic energy for time t
    if n > 1
        for i=1:n
            xnew(i) = 2*x(i) - xold(i) + fx(i)*dt^2;
            ynew(i) = 2*y(i) - yold(i) + fy(i)*dt^2;
            vx(i) = (xnew(i) - xold(i))/(2*dt);
            vy(i) = (ynew(i) - yold(i))/(2*dt);
            k = k + vx(i)^2 + vy(i)^2 ;
        end
    else
        xnew = 2*x - xold + fx*dt^2;
        ynew = 2*y - yold + fy*dt^2;
        vx = (xnew - xold)/(2*dt);
        vy = (ynew - yold)/(2*dt);
        k = k + vx^2 + vy^2 ;
    end
    k = .5*k;
    
    %  create time series of values
    e = k + phi;
    phin(j) = phi/n;
    kn(j) = k/n;
    en(j) = e/n;
    xn(j) = xnew;
    yn(j) = ynew;
    
    % reset positions for next time step
    if n > 1
        for i=1:n
            xold(i) = x(i);
            yold(i) = y(i);
            x(i) = xnew(i);
            y(i) = ynew(i);
        end
    else
        xold = x;
        yold = y;
        x = xnew;
        y = ynew;
    end
    
    % calculate force and energy at new positions for next cycle
    [phi]= phisurf(acon,x,y);
    [fx,fy]= fsurf(acon,x,y);
end