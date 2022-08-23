%%
load(fullfile('..','Results',[mfilename, '_MemoryReentryVaryingFrac.mat']));
figure
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
t = tiledlayout(1,1,'Padding','loose');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 8, 6.8];
ax = nexttile;
plot(memToGCFrac(:,1:end-1), memnum(:,1:end-1), '-o', 'MarkerSize', 4);
xlabel({'Fraction of pre-existing memory', 'cells added to Naive B cell pool'}, 'fontsize', 10)
ylabel('Memory cells', 'fontsize', 10)
leg = legend(ax, {'Dominant', 'Subdominant'}, 'Location', 'east', 'fontsize', 10);
leg.ItemTokenSize = [10,5];
% ylim([10^4, 10^5])
set(ax, 'YScale', 'log')
ylim([10^4, 10^6])
yticks(logspace(3,6,4));

f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_MemoryReentryVaryingFrac.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_MemoryReentryVaryingFrac.pdf']),'ContentType','vector',...
            'BackgroundColor','none');