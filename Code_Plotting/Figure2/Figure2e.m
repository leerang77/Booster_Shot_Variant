%%
load(fullfile('..','Results',[mfilename, '_MeMPCnum.mat']));
figure
t = tiledlayout(1,2,'Padding','tight');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 5.5, 4.7];

ax(1) = nexttile;
egcColors = {[0.3010 0.7450 0.9330], [0.9290 0.6940 0.1250]};
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};

for ep=1:2
    p1(ep) = semilogy(time, memnum{ep}(1,:), 'color', egcColors{ep}, 'LineStyle', '-', 'LineWidth', 1.5);
    hold on
    p2(ep) = semilogy(time(tidx1), memnum{ep}(2,tidx1), 'color', colors{ep}, 'LineStyle', ':', 'LineWidth', 1.5);
    p3(ep) = semilogy(time(tidx2), memnum{ep}(2,tidx2), 'color', colors{ep}, 'LineStyle', '-', 'LineWidth', 1.5);
end

set(gca,'fontname','arial')
xlabel('Time (day)','fontsize',8)
ylabel('Memory cells','fontsize',8)
xticks([28,150])
yticks(logspace(0,6,4))
ylim([1, 10^6])
xlim([0,150])

set(ax(1),'YMinorTick','on')
ax(2) = nexttile;
for ep=1:2
    p1(ep) = semilogy(time, pcnum{ep}(1,:), 'color', egcColors{ep}, 'LineStyle', '-', 'LineWidth', 1.5);
    hold on
    p2(ep) = semilogy(time(tidx1), pcnum{ep}(2,tidx1), 'color', colors{ep}, 'LineWidth', 1.5, 'LineStyle', ':');
    p3(ep) = semilogy(time(tidx2), pcnum{ep}(2,tidx2), 'color', colors{ep}, 'LineStyle', '-', 'LineWidth', 1.5);
    set(gca,'fontname','arial')
    xlabel('Time(day)','fontsize',8)
    ylabel('Plasma cells','fontsize',8)
end
xlim([0,150])
ylim([1, 10^6])
xticks([28,150])
yticks(logspace(0,6,4))
leg = legend([p2(1),p1(1), p3(1),  p2(2), p1(2), p3(2)], {'Vax1 GC Dominant',  'Vax2 EGC Dominant', 'Vax2 GC Dominant','Vax1 GC Subdominant',  ...
    'Vax2 EGC Subdominant', 'Vax2 GC Subdominant'}, 'fontsize', 6, 'NumColumns',2, 'location', 'northoutside');
leg.ItemTokenSize = [10,3];
leg.Layout.Tile = 'North';
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_MeMPCnum.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_MeMPCnum.pdf']),'ContentType','vector',...
            'BackgroundColor','none');