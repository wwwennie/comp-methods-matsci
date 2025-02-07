%% initialize
function[xs,ys]= initsurf2(n,a,nc)
%  pick x and y so they sit near a well
% pick random positions in random wells, so that not too many atoms are
% close to each other
xs = [];
ys = [];

    % place randomly near well min- even random numbers, and then randomly
    %in a repeated cell
%     randx = randi([-a+1,a-1]);
%     randy = randi([-a+1,a-1]);
%     
%     while mod(randx,2*(a/nc)) ~= 0 
%         randx = randi([-a+1,a-1]);
%     end
%     while mod(randy,2*(a/nc)) ~= 0
%         randy = randi([-a+1,a-1]);
%     end
%   
%     x = -0.75+0.5*rand+randx;
%     y = 0.25+0.5*rand+randy;   
% 

% make list of minimum in range of a
% assumes a = 1, form of potential is phisurf
% explore every half coordinate- this is where extrema occur

minsx = [];
minsy = [];

% %make list of minima
% for i = -2*a:2*a
%     xext = 0.5*i; %x-coordinate of extremum
%     for j = -2*a:2*a
%         yext = 0.5*j; % y-coordinate of extremum
% 
%         if abs(fsurf(1,xext,yext)) < 1e-12 % found minimum, with numerical error
%             minsx = [minsx xext];
%             minsy = [minsy yext];
%         end
% 
%     end
% 
% end
  
% minsx
% minsy
% 
% for atom = 1:n
%     randx = randi(length(minsx));
%     randy = randi(length(minsy));
%     
%     % displacement from minimum slightly
%     x = minsx(randx);
%     y = minsy(randy);
%     
%     % append to list of atomic coordinates
%     xs = [xs x];
%     ys = [ys y];
% end

% method 4: brute force 

% make explicit list for nc = 6; half of minima
xnmins = 0.5*[-1 -1 -1 -1 -1 -1 ...
              -5 -5 -5 -5 -5 -5 ...
              -9 -9 -9 -9 -9 -9 ...
              -3 -3 -3 -3 -3 -3 ...
              -7 -7 -7 -7 -7 -7 ...
              -11 -11 -11 -11 -11 -11];
ynmins = 0.5*[1   5  9 -3 -7 -11 ...
              1   5  9 -3 -7 -11 ...
              1   5  9 -3 -7 -11 ...
              3   7 11 -1 -5 -9 ...
              3   7 11 -1 -5 -9 ...
              3   7 11 -1 -5 -9];

xmins = [xnmins -xnmins];
ymins = [ynmins -ynmins];

chosen = [0];
for atom = 1:n
    randind = randi(length(xmins));
    % pick unique sites
    while ~isempty(intersect(chosen,randind))
        randind = randi(length(xmins));
    end
    chosen = [chosen randind];
    % displacement from minimum slightly
    % if plus or minus
    pmx = rand;
    pmy = rand;
    
    if pmx > 0.5
       x = xmins(randind)+0.45*rand;
    else
        x = xmins(randind)-0.45*rand;
    end
    
    if pmy > 0.5
        y = ymins(randind)+0.45*rand;
    else
        y = ymins(randind)-0.45*rand;
    end
    % append to list of atomic coordinates
    xs = [xs x];
    ys = [ys y];
end

%xs=-.75+1.5*rand(n,1);
%ys=-.75+1.5*rand(n,1);
