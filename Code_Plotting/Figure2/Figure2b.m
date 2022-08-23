%%
addpath('..');
load(fullfile('..','Results',[mfilename, '_GCOutcomes.mat']));
figure
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
t = tiledlayout(1,2,'Padding','compact');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 10, 4.5];
subTitles = {'Vax1', 'Vax2'};
for i=1:2
    ax(i) = nexttile;
    p = plotGCNum(param{i}, result{i}, i, nn, colors);
    xlabel(['Time after ', subTitles{i}, ' (day)'], 'fontsize', 8)
    ylim([10^0, 10^4])
    title(subTitles{i}, 'fontsize', 8)
end
leg = legend(ax(1), p(1:2), {'Dominant', 'Subdominant'}, 'Location', 'Northwest', 'fontsize', 6);
leg.ItemTokenSize = [10,5];
ylim([10^0, 10^4])
xticks(ax(2), [0, 60, 120])
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_GCOutcomes.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_GCOutcomes.pdf']),'ContentType','vector',...
            'BackgroundColor','none');