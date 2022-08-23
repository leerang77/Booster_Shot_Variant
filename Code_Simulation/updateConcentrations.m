function [agconc, abconc, Ka, Ka_var] = updateConcentrations(...
    agconc, abconc, Ka, Ka_var, plasmaBlasts, plasmaCells, plasmaCellsEGC, t, param)
%% Documentation
% Summary:
% Updates the antigen and antibody concentrations, antibody affinities

% Outputs: 
%   agconc: 1x3 vector; Concentrations of Antigens; 
%                       soluble, IC1, IC2
%   abconc: 3x2 vector; Concentration of antibodies
%           Dim1 - IgM(natural),IgM(immune),IgG; Dim2-Epitopes
%   Ka: 3x2 vector; Antibody binding affinities (Ka)
%   param: parameter struct
%
%     Note:
%     IC-FDC1: Deposited to FDC by Ig that targets dominant epitope
%     IC-FDC2: Deposited to FDC by Ig that targets subdominant epitope
%     However in this study there is no difference between the two

% Inputs: 
%   agconc, abconc, Ka, Ka_var are same as the outputs
%   plasmaBlasts, plasmaCells, plasmaCellsEGC are arrays of plasma cells
%   t: current time
%   param: parameter struct


%% Function definition
% Check if antigen should be given at current time 
if any(t==param.dose_t)     
   agconc(1) = agconc(1) + param.dose(t==param.dose_t);
end

% Get average affinity of the serum, and soluble Ag and Ab concentrations
R = agconc(1); L = sum(abconc(:));
Ka_avg = sum(sum(Ka.*abconc))/L;

% Get rates for deposition of IC onto FDC 
if L>0 && R>0
    IC_soluble = ((R+L+1/Ka_avg)-sqrt((R+L+1/Ka_avg)^2-4*R*L))/2; %Equil.
    r = param.k_deposit*IC_soluble*(Ka.*abconc)/sum(sum(Ka.*abconc));  
        % 3x2 array; rates of IC-FDC deposition by the Ab species
    if ~isreal(IC_soluble) || isnan(IC_soluble) % Check if im or nan
        error('imaginery number or nan for soluble IC concentration') 
    end
else % If no antigen or antibody
    r = zeros(size(Ka));
end

abconc_old = abconc; %Temporarily save current Ab concentrations

% Check if any of the antibody species concentration will go to zero.
% If this is going to happen, then rescale the reaction rates.
decay_Ab = -r - [0;param.d_IgM;param.d_IgG].*abconc;
rescale_idx = (abconc<-decay_Ab*param.dt);
rescale_factor = abconc./(-decay_Ab*param.dt);
if any(rescale_idx(:))
    fprintf('Reaction rates rescaled. Time t = %.2f', t)
    r(rescale_idx) = r(rescale_idx).*rescale_factor(rescale_idx);
    decay_Ab(rescale_idx) = decay_Ab(rescale_idx).*rescale_factor(rescale_idx);
    if any(isnan(decay_Ab))
        error('Rescaled Ab decay rates contain NaN')
    end
end

% Check if the soluble Ag concentration will go to 0. Rescale if needed.
decay_Ag =  - param.d_Ag*agconc(1)- sum(r(:)); % Net consumption of soluble Ag
if agconc(1)<(-decay_Ag*param.dt)
    rescale_factor = agconc(1)/(-decay_Ag*param.dt);
    decay_Ag = decay_Ag*rescale_factor;
    rnew = r*rescale_factor;
    decay_Ab = decay_Ab + r - rnew;
    r = rnew;
end

%Update the antigen & antibody concentrations for injection / consumption
    agconc(1) = agconc(1) + decay_Ag*param.dt + param.F0*exp(param.k*t)*(t<param.T)*param.dt;
    agconc(2:end) = agconc(2:end) + sum(r,1)*param.dt - agconc(2:end)*param.d_IC*param.dt;
    
    abconc = abconc + decay_Ab*param.dt;
    abconc(abs(abconc)<10^-10)=0;
    if any(agconc(:)<0) || any(abconc(:)<0)
       error('Negative concentration')
    end
        
    
%Update the antibody concentration & affinity for production
M_GC = param.M_GC;
Ig_new = zeros(3,param.n_ep);
Ka_new = repmat({Ig_new}, 1, 2); 
Affinity = cell(2,3);
Target = cell(1,3);

for i=1:2 %WT and Variant
    Affinity{i,1} = plasmaBlasts((i+1)*M_GC+1:(i+2)*M_GC,:);
    Affinity{i,2} = plasmaCells((i+1)*M_GC+1:(i+2)*M_GC,:);
    Affinity{i,3} = plasmaCellsEGC(i+2,:);
    Target{1} = plasmaBlasts(M_GC+1:2*M_GC,:) .* (plasmaBlasts(4*M_GC+1:5*M_GC,:)<(t-param.delay));
    Target{2} = plasmaCells(M_GC+1:2*M_GC,:) .* (plasmaCells(5*M_GC+1:6*M_GC,:)<(t-param.delay));
    Target{3} = plasmaCellsEGC(2,:);
end

for Ig_type = 1:3 %IgM, IgG-GCPC, IgG-EGCPC
    for target=1:param.n_ep
        if any(Target{Ig_type}(:)==target)
            % Ig production
            Ig_new(Ig_type,target) = sum(sum(Target{Ig_type}==target)) * ...
                                (param.r_IgM*(Ig_type==1)+param.r_IgG*(Ig_type>1));
            % Ka of new Ig
            for variant=1:2
               Ka_new{variant}(Ig_type,target) = ...
                   mean(10.^(Affinity{variant,Ig_type}(Target{Ig_type}==target)-9));
            end
        end
    end
end

%Update amounts
abconc(2:3, :) = abconc(2:3, :) + [Ig_new(1,:); sum(Ig_new(2:3, :))]*param.dt;

%update Ka
Ka = {Ka, Ka_var};
for variant=1:2
    current_sum = (abconc_old + param.dt*decay_Ab).*Ka{variant};
    new_sum = current_sum + [0,0,0;1,0,0;0,1,1]*(Ig_new.*Ka_new{variant}*param.dt);
    new_Ka = new_sum./abconc;
    new_Ka(abconc==0)=0;
    if any(new_Ka(:)<0) || any(abs(new_Ka(:))>10^11)
        warning('Error in Ka value: negative or too large')
    end
    new_Ka(isnan(new_Ka)) = 0;
    Ka{variant} = new_Ka;
end
Ka_var = Ka{2};
Ka = Ka{1};
epsilon = 10^-10;
abconc(abconc<epsilon) = 0;
agconc(agconc<epsilon) = 0;
end