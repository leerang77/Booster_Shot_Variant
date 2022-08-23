folders = {'Code_Simulation', 'Code_Parameter_Generation'};
for i=1:length(folders)
addpath(fullfile('..','..',folders{i}))
end

p = base_case_parameters();
param = initializeParameters(p.vaxnum, p.T, p.k, p.numshot, p.E1h, p.dE12, p.p, ...
 p.masking, p.C0, p.w1, p.w2, p.steric, p.memToGCFrac,...
 p.outputprob, p.outputpcfrac, p.rho, p.earlybooster, p.tmax, p.first, p.last);

%%
M_GC = param.M_GC;
E1h = param.E1h;
dE12 = param.dE12;
p = param.p;
f0 = param.f0;
n_mut = param.n_mut;

dE = 0.2;      % class size
classnum = 11;  % total number of class bins
N = 2000;       % total number of naive precursors
NaiveMax = param.NaiveMax;
MemoryMax = param.MemoryReentryMax;

%% Find out the number of B cells in each class
maxclasses = round(linspace(E1h-f0, E1h-f0-dE12, param.n_ep)/dE)+1;% max occupied bin
fitnessarray = linspace(f0, f0+2, classnum); %from 6 to 8 at the interval of 0.2
options = optimoptions('fsolve','Display','off');
r = zeros(param.n_ep,1);
r(1) = fsolve(@(x) (N)-(x^maxclasses(1)-1)/(x-1), 1.1, options);
for i=2:param.n_ep
    if maxclasses(i)>1
        r(i) = fsolve(@(x) (N)-(x^maxclasses(i)-1)/(x-1), 1.1, options);
    else
        r(i) = N;
    end
end
naivebcellsarr = {zeros(param.n_ep,classnum),zeros(param.n_ep,classnum)};
% naivebcellsarr   : n_epx41 array, number of naive B cells in each fitness class
naivebcellsarr{1}(1,1:classnum) = r(1).^(maxclasses(1)-(1:(classnum)));
naivebcellsarr{2}(1,1:classnum) = (1-p*(param.n_ep-1))*r(1).^(maxclasses(1)-(1:(classnum)));
for i=2:param.n_ep
    if maxclasses(i) > 1
        naivebcellsarr{1}(i,1:classnum) = r(i).^(maxclasses(i)-(1:(classnum)));
        naivebcellsarr{2}(i,1:classnum) = p*r(i).^(maxclasses(i)-(1:(classnum)));
    elseif maxclasses(i)==1
        naivebcellsarr{1}(i,1) = N;
        naivebcellsarr{2}(i,1) = p*N;
    end
end
%%
% steric = 0;
% shift = round(-log10(steric)/0.2);
% num(2,1+shift:end) = num(2, 1:end-shift);
% num(2,1:shift) = 0;
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
figure
t = tiledlayout(1,2,'Padding','compact');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 12, 5];
for i=1:2
nexttile;
b = bar(fitnessarray, naivebcellsarr{i}, 'FaceColor', 'flat', 'EdgeColor', 'flat');
for k = 1:2
    b(k).CData = colors{k};
end
set(gca, 'YScale', 'log')
xlabel('Germline affinity (-log_{10}K_d)','fontsize',8)
ylabel('Naive B cell frequency','fontsize',8)
ylim([10^-3.5, 10^3.5])
xticks(6:0.4:8)
yticks(logspace(-3,3,7))
leg = legend({'Dominant', 'Subdominant'}, 'fontsize', 8);
leg.ItemTokenSize = [10,5];
end

% legend({'Ep1', 'Ep2', 'Ep3', 'Ep4'})

f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'.pdf']),'ContentType','vector',...
            'BackgroundColor','none');