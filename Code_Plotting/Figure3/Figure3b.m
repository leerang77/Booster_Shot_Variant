%%
load(fullfile('..','Results',[mfilename, '_AntibodyTiterComparison.mat']));
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
f = figure;
% set(gcf,'Position',[100,100,1000,500])
t = tiledlayout(1,1,'Padding','tight');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 3.5, 4.5];
a1 = nexttile;
x_indx = {[1,2,3]};
n=1;
for i=1:n
% p(1) = plot(x_indx{i},[vax2.WT_titer(1,i),vax2long.WT_titer(1,i),vax3.WT_titer(1,i)], '-o', 'color', 'k', 'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', colors{1});
% hold on
% p(2) = plot(x_indx{i},[vax2.WT_titer(2,i),vax2long.WT_titer(2,i),vax3.WT_titer(2,i)], '-o', 'color', 'k', 'MarkerFaceColor', colors{2}, 'MarkerEdgeColor', colors{2});
p(1) = plot(x_indx{i},[vax2.Var_titer(1,i),vax2long.Var_titer(1,i),vax3.Var_titer(1,i)], '-o', 'color', 'k', 'MarkerSize', 4,'MarkerFaceColor', [1,1,1], 'MarkerEdgeColor', colors{1});
hold on
p(2) = plot(x_indx{i},[vax2.Var_titer(2,i),vax2long.Var_titer(2,i),vax3.Var_titer(2,i)], '-o', 'color', 'k', 'MarkerSize', 4, 'MarkerFaceColor', [1,1,1], 'MarkerEdgeColor', colors{2});
end
set(gca,'YScale','log')
leg = legend(p(1:2), {sprintf('Variant\nDominant'), sprintf('Variant\nSubdominant')},'location','north','fontsize', 6);
leg.ItemTokenSize = [10,3];
% leg.Layout.Tile = 'East';
xlim([0.5,3.5])
ylim(10.^([-3,5]))
yticks(logspace(-3,5,5))
ylabel('Titer (C_{Ab}K_{a})','fontsize',8)
xticks([1,2,3])
xticklabels({'Vax2 1.3m', 'Vax2 5m', 'Vax3 1m'})
set(gca, 'xticklabel', get(a1, 'xticklabel'), 'fontsize', 8)
savefig(f,fullfile('..','figures','Figure3b_Vax3AbTiterDetail.fig') )
exportgraphics(f,fullfile('..','figures','Figure3b_Vax3AbTiterDetail.pdf'),'ContentType','vector',...
            'BackgroundColor','none')


%%
f = figure;
% set(gcf,'Position',[100,100,1000,500])
t = tiledlayout(1,1,'Padding','tight');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 3.5, 4.5];
a1 = nexttile;
x_indx = {[1,2,3]};
for i=1:n
p(1) = plot(x_indx{i},[sum(vax2.WT_titer(:,i)),sum(vax2long.WT_titer(:,i)),sum(vax3.WT_titer(:,i))], '-o', 'MarkerSize', 4, 'color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
hold on
p(2) = plot(x_indx{i},[sum(vax2.Var_titer(:,i)),sum(vax2long.Var_titer(:,i)),sum(vax3.Var_titer(:,i))], '-o', 'MarkerSize', 4, 'color', 'k', 'MarkerFaceColor', [1,1,1], 'MarkerEdgeColor', 'k');
end
set(gca,'YScale','log')
leg = legend(p(1:2), {'WT', 'Variant'},'fontsize', 6, 'location', 'southeast');
leg.ItemTokenSize = [10,3];
% leg.Layout.Tile = 'East';
xlim([0.5,3.5])
ylim(10.^([-3,5]))
yticks(logspace(-3,5,5))
ylabel('Titer (C_{Ab}K_{a})','fontsize',8)
xticks([1,2,3])
xticklabels({'Vax2 1.3m', 'Vax2 5m', 'Vax3 1m'})
set(gca, 'xticklabel', get(a1, 'xticklabel'), 'fontsize', 8)
savefig(f, fullfile('..','figures','Figure3b_Vax3AbTiter.fig'))
exportgraphics(f,fullfile('..','figures','Figure3b_Vax3AbTiter.pdf'),'ContentType','vector',...
            'BackgroundColor','none')