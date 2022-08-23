load(fullfile('..','Results',[mfilename, '_MemoryReentryAffinityNumber.mat']));
figure
t = tiledlayout(1,1,'Padding','tight');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 8, 6];
ax = nexttile;
plot(ctrs, vals, '-*', 'Color', 'k')
hold on
scatter(maxMemAff, memMemoryFraction, 1, [0 0.4470 0.7410], 'filled')
xlabel('Highest affinity of reactivated memory cell (-log_{10}K_d)', 'fontsize', 10)
ylabel({'Fraction of Vax2 GC memory cells', 'from reactivated memory cells'}, 'fontsize', 10)


xlim([6,11])
yticks(0:0.2:1)
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_MemoryReentryAffinityNumber.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_MemoryReentryAffinityNumber.pdf']),'ContentType','vector',...
            'BackgroundColor','none');