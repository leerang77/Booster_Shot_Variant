function agconc_Epmask = epitopeMasking(agconc, abconc, Ka, param)
%% Documentation
% Calculates free antigen concentrations considering epitope masking (or lack
% thereof) based on the current Ag and Ab concentration, Ab affinity

% Output: 
%   agconc_Epmask: 2x2 array; Antigen concentration after epitope masking 
%                  Dim1: soluble, IC-FDC
%                  Dim2: Dominant, Subdominant epitope

% Inputs: 
%   agconc: 1x3 vector; Concentrations of Antigens; 
%                       soluble, IC1, IC2
%   abconc: 3x2 vector; Concentration of antibodies
%           Dim1 - IgM(natural),IgM(immune),IgG; Dim2-Epitopes
%   Ka: 3x2 vector; Antibody binding affinities (Ka)
%   param: parameter struct

%%
masking = param.masking;
agconc_Epmask = zeros(2,param.n_ep); %Dim1: soluble, IC-FDC
                                     %Dim2: Dominant, Subdominant epitope
if sum(agconc)==0
    return
end

if param.n_ep==2
    abconc_steric = sum(abconc) + flip(sum(abconc))*param.steric;
    Ka_avg_steric = (sum(abconc.*Ka) + flip(sum(abconc.*Ka))*param.steric)./abconc_steric;
    abconc = abconc_steric;
    Ka_avg = Ka_avg_steric;
else
    Ka_avg = sum(abconc.*Ka)./sum(abconc); % Mean Ka for each epitope
    abconc = sum(abconc);
end
L = abconc/5; % Ab amount for each epitope
R = sum(agconc); % Total Ag concentration

% Total amount of antigen
IC = (R+L+1./Ka_avg-sqrt((R+L+1./Ka_avg).^2-4*R*L))/2; %covered amounts
IC(isnan(IC)) = 0;

% if param.n_ep == 2
%     x = IC/sum(agconc);
%     IC_steric = sum(agconc)*[x(1)+(1-x(1))*x(2)*param.steric, x(2)+(1-x(2))*x(1)*param.steric];
% else
%     IC_steric = IC;
% end

agconc_Epmask(1,:) = (sum(agconc)-masking*IC)*(agconc(1)/sum(agconc)); %soluble
agconc_Epmask(2,:) = (sum(agconc)-masking*IC)*(sum(agconc(2:end)))/sum(agconc); %FDC
end