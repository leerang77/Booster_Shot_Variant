function numTcells = getnumTcells(param)
%% Documentation
% Summary:
%  Returns the non-dimensionalized availability of helper T cells 
%  over time in a GC

% Outputs:
%  numTcells: array containing number of T cells over time

% Inputs:
%  param: parameters struct
%%
tmax = param.tmax;
d_Tfh = 0.01;
tspan = 0:param.dt:tmax;
numTcells = zeros(size(tspan));
numTmax = param.numTmax;
if tmax<=14
numTcells = numTmax*(tspan/14);
else
d14idx = round(14/param.dt+1);
numTcells(1:d14idx) = numTmax*(tspan(1:d14idx)/14);
for i=d14idx+1:length(tspan)
   numTcells(i) = numTcells(i-1)*exp(-d_Tfh*param.dt); 
end
end
end