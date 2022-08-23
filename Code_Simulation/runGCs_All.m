function runGCs_All(varargin)
%% Documentation
% Summary:
%  Get all parameters as inputs except for vaxnum and tmax. Then,
%  sequentially run the standard Vax 1 through Vax 3 procedure. 

[T, k, numshot, E1h, dE12, p, ...
 masking, C0, w1, w2, steric, memToGCFrac,...
 outputprob, outputpcfrac, rho, earlybooster, first, last] = deal(varargin{:});

fprintf('Vax1')
vaxnum = 1; tmax = 28;
result = runGCs(vaxnum, T, k, numshot, E1h, dE12, p, ...
 masking, C0, w1, w2, steric, memToGCFrac,...
 outputprob, outputpcfrac, rho, earlybooster, tmax, first, last)

clear result

pause(30);

fprintf('Vax2')
vaxnum = 2; tmax = 180;
result = runGCs(vaxnum, T, k, numshot, E1h, dE12, p, ...
 masking, C0, w1, w2, steric, memToGCFrac,...
 outputprob, outputpcfrac, rho, earlybooster, tmax, first, last)

clear result

pause(30);

fprintf('Vax3')
vaxnum = 3; tmax = 28;
result = runGCs(vaxnum, T, k, numshot, E1h, dE12, p, ...
 masking, C0, w1, w2, steric, memToGCFrac,...
 outputprob, outputpcfrac, rho, earlybooster, tmax, first, last)
end