% simple lattice sum for force with cutoffs and 
% minimum image convention, plus calculating radial distribution function
%
%  we calculate force (fx,fy,fz), energy (u), and
%     part of the pressure (w)
function[u,w,fx,fy,fz,bp,ng]= fLJsumg(a,n,x,y,z,nbin)
%  set force components, potential energy, and pressure to 0
fx=zeros(n,1);
fy=zeros(n,1);
fz=zeros(n,1);
u = 0;
w = 0;


%part for computing g(r)
rc = a/2;
xb = rc/nbin;
g = zeros(nbin,1);
bp = zeros(nbin,1);
ng = zeros(nbin,1);

for i = 1:n-1  % note limits
    ftx = 0;
    fty = 0;
    ftz = 0;
    for j=i+1:n  %  note limits
% mimimum image convention
        dx = x(j) - x(i);
        dy = y(j) - y(i); 
        dz = z(j) - z(i);
        dx = dx - round(dx);
        dy = dy - round(dy);
        dz = dz - round(dz);
        dist = a*sqrt(dx^2 + dy^2 + dz^2);
        if dist <= rc
            dphi = (2/dist^(12)-1/dist^6);
            ffx = dphi*a*dx/dist^2;
            ffy = dphi*a*dy/dist^2;
            ffz = dphi*a*dz/dist^2;
            ftx = ftx + ffx;
            fty = fty + ffy;
            ftz = ftz + ffz;
            phi = (1/dist^(12)-1/dist^6);
            u = u +  phi;
            w = w + dphi;
%  add -f to sum of force on j
            fx(j) = fx(j) - ffx;
            fy(j) = fy(j) - ffy;
            fz(j) = fz(j) - ffz;
            
            %part for g(r)
            ib = floor(dist/xb);
            g(ib) = g(ib) + 1;
        end
    end
  %  sum up force on i (fi)
    fx(i) = fx(i) + ftx;
    fy(i) = fy(i) + fty;
    fz(i) = fz(i) + ftz;
end
%  need to multiply LJ by 4 and force and pressure by 24
%  also need to correct sign in f 
u = 4*u;
w = 24*w;
for i=1:n
    fx(i) = -24*fx(i);
    fy(i) = -24*fy(i);
    fz(i) = -24*fz(i);
end

%for computing g(r)
% normalize and create proper distances
factor = 2*a^3/(4*pi*n^2*xb);
for i=1:nbin
    bp(i) = (i+1/2)*xb;
    ng(i) = factor*g(i)/(i*xb)^2;
end