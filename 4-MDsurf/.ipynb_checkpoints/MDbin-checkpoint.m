% MD code for binning
% corresponds to problem 1, part c

% Input: data = contains data to find average
%        nbins = # bins to average over
% Ouput: avg = average of data

function avg = MDbin(data,nbin)

nsteps = length(data);
binsteps = nsteps/nbin;  %width of bin
binavgs = [];

% for ease of calculation, make data divisible by nbin
if mod(nsteps,nbin) ~= 0 
    error('Data size and # bins not divisible by each other!!!')
end

%bin averages
idx = 1;
for i = 1:nbin
   binavg = 0;
   for j = 1:binsteps
        binavg = binavg + data(idx);
        idx = idx + 1;
   end
   
   binavgs = [binavgs 1/binsteps*binavg];
end

%data averages
avg = 1/nbin*sum(binavgs);

end