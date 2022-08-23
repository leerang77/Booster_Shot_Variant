%%
% cd(fileparts(which(mfilename)));
load(fullfile('..','Results',[mfilename, '_AgConc.mat']));
figure
t = tiledlayout(1,1,'Padding','tight');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 4.5, 3.8];
ax = nexttile;
p1 = semilogy(agconc_x, agconc_y(1,:), 'color', [0, 0, 0, 1], 'LineStyle', ':');
hold on
p2 = semilogy(agconc_x, agconc_y(2,:), 'color', [0, 0, 0, 1], 'LineStyle', '-');
set(gca,'fontname','arial')
xlabel('Time (day)','fontsize',8)
ylabel('Antigen concentration (nM)','fontsize',8)
leg = legend(ax, [p1, p2], {'soluble', 'IC-FDC'},'Location','NorthEast','Orientation','Vertical', 'fontsize', 6);
leg.ItemTokenSize = [10,5];
% leg.Layout.Tile = 'North';

xlim([0,30*3])
% xticks(0:30:90)
ylim([5*10^-6, 10])
yticks(logspace(-5,1,7));
labels = cell(1,7);
labels(1:2:7) = repmat({''},1,4);
labels(2:2:6) = cellfun(@(x) sprintf('10^{%d}',x), {-4,-2,0},'UniformOutput',false);
yticklabels(labels);

f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_AgConc.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_AgConc.pdf']),'ContentType','vector',...
            'BackgroundColor','none');