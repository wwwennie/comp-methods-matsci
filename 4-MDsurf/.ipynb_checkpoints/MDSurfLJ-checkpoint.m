% MD code for 2D system of surface atoms using Verlet Algorithm + LJ pot
%  input:  n = number of atoms
%          a = simulation cell size
%          Tin = temperature for initializing velocities
%          acon = parameter in surface potential; strength of interaction
%          nsteps = number of time steps
%          dt = time step
%  output:  un = potential energy at each time step
%           kn = kinetic energy at each time step
%           en = total energy at each time step
%           tn = temperature at each time step
%           pn = pressure at each time step

function[phin,kn,en,xn,yn]= MDSurfLJ(n,a,nc,Tin,acon,rc,nsteps,dt)

% switch to temperature base for initial velocities
% use initsurf for initial atomic positions only
[xs,ys] = initsurf2(n,a,nc);

% use Maxwell-Distribution for velocities
[vxs,vys,pxs,pys] = MBinit(n,Tin);

% calculate initial energy and forces
[phi]= phisurf(acon,xs,ys);
[fxsurf,fysurf]= fsurf(acon,xs,ys);

%need input scaled coordinates
xscaled = xs./a;
yscaled = ys./a;
[ulj,wlj,fxlj,fylj] = fLJ2D(a,n,rc,xscaled,yscaled);

netfx = fxsurf+fxlj;
netfy = fysurf+fylj;
netpot = phi+ulj;

% now start the time stepping with the verlet algorithm
% initialize variables
xold = zeros(n,1);
yold = zeros(n,1);

xnew = zeros(n,1);
ynew = zeros(n,1);

% first find the positions at t-dt
for i=1:n
    xold(i) = xs(i) - vxs(i)*dt + .5*netfx(i)*dt^2;
    yold(i) = ys(i) - vys(i)*dt + .5*netfy(i)*dt^2;
end


% start the time steps
for j=1:nsteps
    k = 0;
 %  find positions for time t + dt
 %  find velocities for time t
 %  find kinetic energy for time t
    
    for i=1:n
        xnew(i) = 2*xs(i) - xold(i) + netfx(i)*dt^2;
        ynew(i) = 2*ys(i) - yold(i) + netfy(i)*dt^2;
        vxs(i) = (xnew(i) - xold(i))/(2*dt);
        vys(i) = (ynew(i) - yold(i))/(2*dt);
        k = k + vxs(i)^2 + vys(i)^2 ;
    end
    
    k = .5*k;
    
    %  create time series of values
    e = k + netpot;
    phin(j,:) = netpot/n;
    kn(j,:) = k/n;
    en(j,:) = e/n;
    xn(j,:) = xnew;
    yn(j,:) = ynew;
    
    % reset positions for next time step
 
    for i=1:n
        xold(i) = xs(i);
        yold(i) = ys(i);
        xs(i) = xnew(i);
        ys(i) = ynew(i);
    end

    
    % calculate force and energy at new positions for next cycle
    [phi]= phisurf(acon,xs,ys);
    xscaled = xs./a;
    yscaled = ys./a;
    [fxsurf,fysurf]= fsurf(acon,xs,ys);
    [ulj,wlj,fxlj,fylj] = fLJ2D(a,n,rc,xscaled,yscaled);
    netfx = fxsurf+fxlj;
    netfy = fysurf+fylj;
    netpot = phi+ulj;
end