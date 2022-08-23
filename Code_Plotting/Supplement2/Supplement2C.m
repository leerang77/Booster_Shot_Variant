%%
load(fullfile('..','Results',[mfilename, '_Unique_Lineages_Size.mat']))
figure
set(gca,'YScale','log')
xlabel('rank')
ylabel('size')
set(gcf,'Position',[100,100,300,300])
t = tiledlayout(1,1,'Padding','tight');
t.Units = 'centimeters';
t.InnerPosition = [3, 3, 3.7, 3.7];
ax = nexttile;
colors = {[117,107,177]/256,[253,174,97]/256,[28,144,153]/256};
plot(1:length(GC{1}), GC{1}, 'Linewidth', 1.5, 'Color', colors{1})
hold on
plot(1:length(EGC{2}), EGC{2}, 'Linewidth', 1.5, 'Color', colors{2});
plot(1:length(GC{2}), GC{2}, 'Linewidth', 1.5, 'Color', colors{3});

colormap(ax,'turbo')
set(ax, 'YScale', 'log')
leg = legend({'Vax1 1m GC', 'Vax2 5m EGC', 'Vax2 5m GC'}, 'Location', 'NorthEast');
leg.ItemTokenSize = [15, 5];
xlabel('Rank of lineage size')
ylabel('Lineage size')
xlim([0,500])
yticks(logspace(0,5,6));
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_Unique_Lineages_Size.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_Unique_Lineages_Size.pdf']),'ContentType','vector',...
            'BackgroundColor','none');