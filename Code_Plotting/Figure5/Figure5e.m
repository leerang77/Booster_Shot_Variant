%%
load(fullfile('..','Results',[mfilename, '_StericChange.mat']))
colors = {[0 0.4470 0.7410],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980]};
colors = {colors{1}, colors{2:param{1}.n_ep-1}, colors{end}};
figure
t = tiledlayout(1,1,'Padding','tight');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 5.05, 4.8];
ax = nexttile;
plot(sterics, memnum, '-o', 'MarkerSize', 4);
xlabel({'Epitope Overlap'}, 'fontsize', 8)
ylabel('Memory cells', 'fontsize', 8)
leg = legend(ax, {'Dominant', 'Subdominant'}, 'Location', 'southwest', 'fontsize', 8);
leg.ItemTokenSize = [10,5];
% ylim([10^4, 10^5])
set(ax, 'YScale', 'log')
xlim([0.2,0.7])
ylim([10^4, 10^6])
yticks(logspace(3,6,4));
xticks(0.2:0.1:0.7)

f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_StericChange.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_StericChange.pdf']),'ContentType','vector',...
            'BackgroundColor','none');