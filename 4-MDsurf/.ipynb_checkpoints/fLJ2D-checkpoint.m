% simple lattice sum for force with cutoffs and 
% minimum image convention for 2D surface potential
%
%  we calculate force (fx,fy,fz), energy (u), and
%     part of the pressure (w)
function[u,w,fx,fy]= fLJ2D(a,n,rc,x,y)
    %  set force components, potential energy, and pressure to 0
    fx=zeros(1,n);
    fy=zeros(1,n);

    u = 0;
    w = 0;
    for i = 1:n-1  % note limits
        ftx = 0;
        fty = 0;

        for j=i+1:n  %  note limits
    % mimimum image convention
            dx = x(j) - x(i);
            dy = y(j) - y(i); 

            dx = dx - round(dx);
            dy = dy - round(dy);

            dist = a*sqrt(dx^2 + dy^2);
            if dist <= rc
                dphi = (2/dist^(12)-1/dist^6);
                ffx = dphi*a*dx/dist^2;
                ffy = dphi*a*dy/dist^2;

                ftx = ftx + ffx;
                fty = fty + ffy;

                phi = 4*(1/dist^(12)-1/dist^6);
                u = u +  phi;
                w = w + dphi;
    %  add -f to sum of force on j
                fx(j) = fx(j) - ffx;
                fy(j) = fy(j) - ffy;

            end
        end
      %  sum up force on i (fi)
        fx(i) = fx(i) + ftx;
        fy(i) = fy(i) + fty;

    end
    %  need to multiply LJ by 4 and force and pressure by 24
    %  also need to correct sign in f 
    u = 4*u;
    w = 24*w;
    for i=1:n
        fx(i) = -24*fx(i);
        fy(i) = -24*fy(i);

    end
end