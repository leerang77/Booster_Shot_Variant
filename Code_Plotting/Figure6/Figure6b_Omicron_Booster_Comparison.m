fprintf(mfilename)
M_GC = 2000;
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
addpath(fullfile('..','Code_Parameter_Generation'));
addpath(fullfile('..','Code_Result_Analysis'));
% Load prime result, no Ep masking


result = cell(5,2); %Row: WT booster, Omicron Booster
                    %Column: Ep masking 0.3, 0.7; Rho 0.9, p=0.5; memToGC0.04; 1; 
                    
tic
%% Ep masking 0.3
[T, k, numshot, E1h, dE12, p, n_ep, ...
    masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
    w1, w2, tmax, first, last, numfrag, steric, memToGC] = base_case_parameters();
vaxnum = 3;
tmax = 28;
masking = 1;
steric = 0.3;
result{1,1} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
    masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
    w1, w2, tmax, first, last, steric, memToGC);
result{1,2} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
    masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
    w1, w2, tmax, first, last, steric, memToGC, 'Omicron_Booster');
toc
%% Ep masking 0.7
steric = 0.7;
result{2,1} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
    masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
    w1, w2, tmax, first, last, steric, memToGC);

result{2,2} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
    masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
    w1, w2, tmax, first, last, steric, memToGC, 'Omicron_Booster');
toc
%% Rho=0.9
masking = 0;
steric = 0;
result{3,1} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
    masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
    w1, w2, tmax, first, last, steric, memToGC);

result{3,2} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
    masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
    w1, w2, tmax, first, last, steric, memToGC, 'Omicron_Booster');
toc
%% E1h=7.4,dE12=0.4,p=0.5
E1h = 7.4;
dE12 = 0.4;
p = 0.5;
last = 1000;
result{4,1} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
    masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
    w1, w2, tmax, first, last, steric, memToGC);

result{4,2} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
    masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
    w1, w2, tmax, first, last, steric, memToGC, 'Omicron_Booster');
toc
%% memTOGC fraction 0.04

% memToGC = 0.04;
% result{5,1} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
%     masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
%     w1, w2, tmax, first, last, steric, memToGC);
% % 
% result{5,2} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
%     masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
%     w1, w2, tmax, first, last, steric, memToGC, 'Omicron_Booster');
% toc

%% rho 1
E1h = 7;
dE12 = 0.4;
p = 0.1;
last = 2000;
memToGC = 0;
result{5,1} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
    masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
    w1, w2, tmax, first, last, steric, memToGC, 'Rho1');

result{5,2} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
    masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
    w1, w2, tmax, first, last, steric, memToGC, fullfile('Rho1','Omicron_Booster'));
toc
%%
summary = cell(size(result));
[n,m] = size(result);
WT_titer = zeros(2,n,m);
Var_titer = zeros(2,n,m);
for i=1:n
    for j=1:m
        summary{i,j} = AbConcentration(result{i,j}, result{i,j}.param);
        WT_titer(:,i,j) = (summary{i,j}.titer_geomean{1}(:,end));
        Var_titer(:,i,j) = (summary{i,j}.titer_geomean{2}(:,end));
    end
end
x=1:n;
save('Figure6b_Omicron_Booster_Comparison.mat', 'summary', 'WT_titer', 'Var_titer');


%% Figure
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
load('Figure6b_Omicron_Booster_Comparison.mat')

orders = [3,4,2,1,5];
f=figure;
t = tiledlayout(1,1,'Padding','tight');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 9, 6];
% for i=1:length(pairs)
% % plot(i*2-1:i*2, sum(squeeze(WT_titer(:,pairs{i}(1),:))), '-o', 'color', 'k', 'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', colors{1}) % 0.3, WT booster and omicron booster
% hold on
% plot(i*2-1:i*2, sum(squeeze(Var_titer(:,pairs{i}(1),:))), '-o', 'color', 'k','MarkerFaceColor', colors{2}, 'MarkerEdgeColor', colors{2}) % 0.3, WT booster and omicron booster
% % plot(i*2-1:i*2, sum(squeeze(WT_titer(:,pairs{i}(2),:))), '-o', 'color', 'k','MarkerFaceColor', [1,1,1], 'MarkerEdgeColor', colors{1}) % 0.3, WT booster and omicron booster
% hold on
% plot(i*2-1:i*2, sum(squeeze(Var_titer(:,pairs{i}(2),:))), '-o', 'color', 'k','MarkerFaceColor', [1,1,1], 'MarkerEdgeColor', colors{2}) % 0.3, WT booster and omicron booster
% end

vals = zeros(5,2);
for i=1:length(orders)
    vals(i,:) = [sum(squeeze(Var_titer(:,orders(i),:)))];
end
nexttile;
% plot(i*2-1:i*2, sum(squeeze(WT_titer(:,pairs{i}(1),:))), '-o', 'color', 'k', 'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', colors{1}) % 0.3, WT booster and omicron booster
b = bar(1:5, vals, 'FaceColor', 'flat', 'EdgeColor', 'r') % 0.3, WT booster and omicron booster
b(1).CData = colors{2};
b(2).CData = [1,1,1];
% plot(i*2-1:i*2, sum(squeeze(WT_titer(:,pairs{i}(2),:))), '-o', 'color', 'k','MarkerFaceColor', [1,1,1], 'MarkerEdgeColor', colors{1}) % 0.3, WT booster and omicron booster
% bar(i*2-1:i*2, sum(squeeze(Var_titer(:,pairs{i}(2),:)))) % 0.3, WT booster and omicron booster
set(gca,'Yscale','log')
ylim([10^0, 10^5])
yticks(logspace(0,4,5))
xlim([0, 6])
xticks([1,2,3,4,5])
xticklabels({'Base', 'Precursor50', 'Masking70',...
    'Masking30', 'Conservation100'})
vals(:,2)./vals(:,1)
ylabel('1m Vax3 variant titer', 'fontsize', 8)
legend({'WT Booster', 'Variant Booster'}, 'Location', 'NorthWest', 'fontsize', 8)
set(gca,'Fontsize',8)
savefig(f,fullfile('figures','Figure6b_VaxWithVariant.fig') )
exportgraphics(f,fullfile('figures','Figure6b_VaxWithVariant.pdf'),'ContentType','vector',...
            'BackgroundColor','none')
