%%
load(fullfile('..','Results',[mfilename, '_MemoryByInitialAffinity.mat']));
figure
t = tiledlayout(1,2,'Padding','tight');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 4.6, 4.84];
subplotnames = {'Vax1 1m', 'Vax2 5m'};
num_sim = 10;
for i=1:2
    ax(i) = nexttile;
    plotbar(num_by_E0_sum{i}/num_sim, ax(i))
    set(gca,'fontname','arial')
    set(gca, 'FontSize', 8)
    title({subplotnames{i};''})     
end
ylim(ax(1),[0,3*10^3])
ylim(ax(2),[0,2*10^5])
% xlim(ax(1),[5.8,10])
% xlim(ax(2),[5.8,10])
ax(1).YAxis.Exponent = 3;
ax(2).YAxis.Exponent = 5;
xlabel(t,{'Germline affinity (-log_{10}K_{d})'}, 'FontSize', 8)
ylabel(t,'Memory cells', 'FontSize', 8)
leg = legend(ax(2), {'Dominant', 'Subdominant'}, 'location','NorthOutside','Orientation','Horizontal', 'fontsize', 6);
leg.ItemTokenSize = [10,3];
leg.Layout.Tile = 'North';
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_MemoryByInitialAffinity.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_MemoryByInitialAffinity.pdf']),'ContentType','vector',...
            'BackgroundColor','none');
for i=1:2
    num_by_E0_sum{i}(isnan(num_by_E0_sum{i}))=0;
    sum(sum(num_by_E0_sum{i}(:,6:end)))/sum(sum(num_by_E0_sum{i}))
end

function b = plotbar(num_by_E0_sum, ax)
barcolors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.3010 0.7450 0.9330], [0.9290 0.6940 0.1250]};
b = bar(ax, (6:0.2:8), num_by_E0_sum', 'stacked', 'FaceColor', 'flat', 'EdgeColor', 'flat', 'barwidth', 0.7);
b(1).CData = barcolors{1};
b(1).EdgeColor = barcolors{1};
b(2).CData = barcolors{2};
b(2).EdgeColor = barcolors{2};
end