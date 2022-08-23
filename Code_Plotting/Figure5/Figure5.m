Figure5a_Epitope_Masking_Visualization
Figure5b_Epitope_Blocking_Data
%%
clear all
addpath(fullfile('..','Code_Parameter_Generation'));
addpath(fullfile('..','Code_Result_Analysis'));
% Load prime result, no Ep masking
[T, k, numshot, E1h, dE12, p, n_ep, ...
    masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
    w1, w2, tmax, first, last, numfrag, steric, memToGC] = base_case_parameters();
result = cell(1,2);
% vaxnum = 3;
% tmax = 28;
% steric = 0;
% memToGC = 0;
% masking = 1;
% result{1} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
%     masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
%     w1, w2, tmax, first, last, steric, memToGC);
vaxnum = 2;
tmax = 180;
steric = 0.3;
masking = 1;
result{2} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
    masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
    w1, w2, tmax, first, last, steric, memToGC, 'Data');
% param{1} = result{1}.param;
param{2} = result{2}.param;
M_GC = 2000;
colors = {[0 0.4470 0.7410],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980]};
colors = {colors{1}, colors{2:param{2}.n_ep-1}, colors{end}};

%%
Figure5c_Epitope_Masking_Steric_Hindrance_Agconc
Figure5d_Epitope_Masking_Steric_Hindrance_GC

%%
[T, k, numshot, E1h, dE12, p, n_ep, ...
    masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
    w1, w2, tmax, first, last, numfrag, steric, memToGC] = base_case_parameters();
vaxnum = 3;
tmax = 28;
steric = 0;
memToGC = 0;
result{1} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
    masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
    w1, w2, tmax, first, last, steric, memToGC,'Data');
steric = 0.3;
masking = 1;
result{2} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
    masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
    w1, w2, tmax, first, last, steric, memToGC,'Data');

Figure5e_NoMaskingvsMasking
