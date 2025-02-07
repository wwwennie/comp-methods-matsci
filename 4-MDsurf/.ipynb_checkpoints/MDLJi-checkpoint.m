% MD code for 3D system of LJ atoms with initialization modification
%  input:  nc = number of cells of fcc unit cell
%          density = density of LJ system;
%          tin = initial temperature    
%          nsteps = number of time steps
%          dt = time step
 
%  output:  un = potential energy at each time step
%           kn = kinetic energy at each time step
%           en = total energy at each time step
%           tn = temperature at each time step
%           pn = pressure at each time step


function[un,kn,en,tn,pn,a,x,y,z,vx,vy,vz]= MDLJi(density,nsteps,dt,n,x,y,z,vx,vy,vz)

% calculate some useful quantities
vol = n/density;
a = vol^(1/3);
rc = a/2;

% calculate initial energy and forces
[u,w,fx,fy,fz] = fLJsum(a,n,rc,x,y,z);

% now start the time stepping with the verlet algorithm
% initialize variables
xold = zeros(n,1);
yold = zeros(n,1);
zold = zeros(n,1);
xnew = zeros(n,1);
ynew = zeros(n,1);
znew = zeros(n,1);

% first find the positions at t-dt
for i=1:n
    xold(i) = x(i) - vx(i)*dt/a + .5*fx(i)*dt^2/a;
    yold(i) = y(i) - vy(i)*dt/a + .5*fy(i)*dt^2/a;
    zold(i) = z(i) - vz(i)*dt/a + .5*fz(i)*dt^2/a;
end

% start the time steps
for j=1:nsteps
    k = 0;
 %  find positions for time t + dt 
 %  find velocities for time t
 %  find kinetic energy for time t
    for i=1:n
        xnew(i) = 2*x(i) - xold(i) + fx(i)*dt^2/a;
        ynew(i) = 2*y(i) - yold(i) + fy(i)*dt^2/a;
        znew(i) = 2*z(i) - zold(i) + fz(i)*dt^2/a;
        vx(i) = a*(xnew(i) - xold(i))/(2*dt);
        vy(i) = a*(ynew(i) - yold(i))/(2*dt);
        vz(i) = a*(znew(i) - zold(i))/(2*dt);
        k = k + vx(i)^2 + vy(i)^2 + vz(i)^2;
    end
    k = .5*k;
    temp = 2*k/(3*n);
%  create time series of values
    e = k + u;
    un(j) = u/n;
    kn(j) = k/n;
    en(j) = e/n;
    tn(j) = temp;
    pn(j) = density*temp + w/(3*vol);
        
% reset positions for next time step
    for i=1:n
        xold(i) = x(i);
        yold(i) = y(i);
        zold(i) = z(i);
        x(i) = xnew(i);
        y(i) = ynew(i);
        z(i) = znew(i);
    end
% calculate force and energy at new positions for next cycle
    [u,w,fx,fy,fz] = fLJsum(a,n,rc,x,y,z);
end
