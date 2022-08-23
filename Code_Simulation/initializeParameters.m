function param = initializeParameters(varargin)
%% Documentation
% Summary:
%   Initializes the variables based on inputs.
%   Also initialize parameters whose values are fixed. 

% Inputs:
%  varargin: cell array of the parameters in order (see lines 11-13)
%   vaxnum: Integer between 1 and 4; indicates dose number
%   T: 0 indicates bolus
%   k: 0 indicates bolus
%   numshot: 1 indicates bolus
%   E1h: Between 6 and 8. Defines naive B cell germline affinities
%   dE12: Between 0 and 1. Defines naive B cell germline affinities
%   p: Between 0 and 1. Fraction of naive B cells that are subdominant
%   masking: 1 if epitope masking, 0 if not
%   C0: Reference antigen concentration
%   w1: Typically 0. If not, then alternative Ag capture model is used.
%       Defines saturation of antigen capture 
%   w2: Selection stringency. Typically 0.5; varied between 0.3 and 1
%   steric: 0 if no masking. Varied between 0 and 1 if masking is 1
%   memToGCFrac: Fraction of pre-existing memory cells that can enter GC
%   outputprob: Fraction of selected GC B cells that exit
%   outputpcfrac: Fraction of exiting GC B cells that become plasma cell
%   rho: Epitope conservation level for subdominant. Typically 0.95
%   earlybooster: Typically 0. Only for 3rd dose. If not 0, then it defines
%   how many days after Vax2 the Vax3 is given (eg. 42)
%   tmax: Time duration of simulation 
%   first: First index of GCs (eg. 201)
%   last: Last index of GCs (eg. 400) 

% Output:
%   param: Struct having parameters and constants as fields
%          See the body of the code for the fields and explanations
%% Initialization
[vaxnum, T, k, numshot, E1h, dE12, p, ...
 masking, C0, w1, w2, steric, memToGCFrac,...
 outputprob, outputpcfrac, rho, earlybooster, tmax, first, last] ...
 = deal(varargin{:});
param = struct();

%% Simulation numbers 
param.first = first; 
param.last = last;
param.M_GC = last-first+1; % Number of GCs;

%% Data storage
param.N_GC_MAX = 3000; % Max B cell number per GC
param.N_PC_MAX = 10000; % Max PC number per GC 
param.N_PB_MAX = 1000; % Max PB number per GC 
% (However PB production is not considered in this simulation) 
param.naivefieldnum = 7; % Number of properties that naive B cells have 
param.gcfieldnum = 5; % Number of properties that GC B cells have
param.pcfieldnum = 6; % Number of properties that PCs have
param.memfieldnum = 9; % Number of properties that memory cells have

%% Dosing scheme
param.vaxnum = vaxnum;
param.T = T;
param.k = k;
param.numshot = numshot;
param.earlybooster = earlybooster;

%% Naive B cell affinity Distribution
param.E1h = E1h;
param.dE12 = dE12;
param.p = p;
param.n_ep = 2; % Number of epitopes. Dominant and subdominant epitope
param.NaiveMax = 2010; % Max number of Naive B cells per GC 

%% Re-entry of memory B cells into GCs
param.memToGCFrac = memToGCFrac;
param.MemoryReentryMax = 200; % Max number of pre-existing memory cells 
                              % that can re-enter GC, per GC

%% Antibody Production and Effect
param.production = 1; % Parameter that is used to change the rate 
                      % of antibody production. Fixed in this simulation
param.delay = 2;
param.masking = masking;
param.steric = steric;
param.outputprob = outputprob;
param.outputpcfrac = outputpcfrac;
param.mutprob = 1; % Probability of mutation of a daughter cell. Fixed 
param.IgM0 = 0.01; % nM; Initial IgM amount 
param.r_IgG = param.IgM0*param.production/param.M_GC; %nM per PC per day;
                           %IgG production rate normalized by number of GCs
param.r_IgM = param.IgM0*param.production/param.M_GC; %nM per PB per day;
                           %IgM production rate normalized by number of GCs

%% Antigen Dose and Effect
param.Ag0 = 10; % nM; Initial soluble Ag amount
param.Ageff = 0.01; % Effectiveness of soluble Ag at B cell activation
param.C0 = C0;
param.F0 = nan; % Not used in this simulation; for slow delivery
param.dose = []; % Not used in this simulation; for slow delivery
param.dose_t = []; % Not used in this simulation; for slow delivery

%% Time
param.tmax = tmax;
param.dt = 0.01; % day; time interval of simulation
param.tspan_summary = 0:(25*param.dt):param.tmax; % time points at which
                               % statistics from the simulation are saved

%% Activation and Competition Parameters
param.f0 = 6; % Reference binding affinity 
param.activation_threshold = 1; % Amount of antigen captured for activation
param.w1 = w1;
param.w2 = w2;

%% Mutation Parameters
param.MutationPDF = [3.1, 1.2, 3.08];
param.rho = rho;
param.n_res = 80; % Number of residues that each B cell has

%% Reaction Rates
param.k_deposit = 24; % per day; Rate at which IC is transported to FDC
param.d_Ag = 3; % per day; Rate of antigen decay 
param.d_IC = 0.15; % per day; Rate of decay for IC on FDC
param.d_IgM = log(2)/28; % per day; Decay rate of IgM; Half-life 28 days 
param.d_IgG = log(2)/28; %per day; Decay rate of IgG; Half-life 28 days
param.d_pc = log(2)/4; %per day; Decay rate of PCs; Half-life 4 days

%% Recruitment, birth, death rates
param.lambdamax = 1; % per day; Max rate of GC entry for naive B cells
param.betamax = 2.5; % per day; Max rate of GC B cell positive selection.
                      % Division of up to ~4 times a day is possible
param.betamaxEGC = 2.5; % per day
param.mu = 0.5; % per day; Death rate of GC B cells 
param.Nmax = 10; % per day; Defines number of naive B cells that enter GC
param.numTmax = 1200; % Max non-dimensionalized helper T cell availability
param.naiveprolnum = 4; % Number of copies made when Naive B cell enters GC
end