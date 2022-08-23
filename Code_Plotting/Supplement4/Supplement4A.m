%% Vax1,2
suffixes = {'Vax1', 'Vax2'};
f = figure();
t = tiledlayout(1,2,'Padding','tight');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 7.5, 3.5];
ax = cell(1,2);
for fignum=1:2
suffix = suffixes{fignum};
load(fullfile('..','Results',[mfilename,'_',suffix,'_VaryingC0.mat']));
ax{fignum} = nexttile;
plot(C0s, memnum, '-o', 'MarkerSize', 4);
set(ax{fignum}, 'XScale', 'log')
set(ax{fignum}, 'YScale', 'log')
yticks(logspace(1,6,6))
if fignum==2
    leg = legend({'Dominant', 'Subdominant'}, 'Location', 'Southeast', 'fontsize', 6);
    leg.ItemTokenSize = [15,5];
end
xlabel({'C_0'}, 'fontsize', 8)
end


ylabel(t,'Memory cells', 'fontsize', 8)

title(ax{1},'Vax1 1m GC')
title(ax{2},'Vax2 5m GC')
ylim(ax{1}, [10, 10^5])
ylim(ax{2}, [100,10^6])
savefig(f, fullfile('..','figures',[mfilename,'_VaryingC0.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_VaryingC0.pdf']),'ContentType','vector',...
            'BackgroundColor','none');