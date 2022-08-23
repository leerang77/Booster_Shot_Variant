%%
load(fullfile('..','Results',[mfilename, '_Antibody.mat']));
figure
t = tiledlayout(1,1,'Padding','compact');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 5, 5];
nexttile;
idx = 1:4:833;
semilogy(time(idx), WT_titer(1,(idx)), 'Color', colors{1}, 'LineWidth', 1.5)
hold on
semilogy(time(idx), WT_titer(2,(idx)), 'Color', colors{2}, 'LineWidth', 1.5)
semilogy(time(idx), Variant_titer(1,(idx)), ':', 'Color', colors{1}, 'LineWidth', 1.5)
semilogy(time(idx), Variant_titer(2,(idx)), ':', 'Color', colors{2}, 'LineWidth', 1.5)
ylim([10^-2, 10^3])
yticks(logspace(-2,4,4))
xticks([0,28,90,180])
leg = legend('WT-Dominant', 'WT-Subdominant', 'Variant-Dominant', 'Variant-Subdominant', 'fontsize', 6, 'NumColumns', 2);
leg.ItemTokenSize = [10,3];
leg.Layout.Tile = 'North';
set(gca,'fontname','arial')
xlabel('Time (day)','fontsize',8)
ylabel('Titer (C_{Ab}K_{a})','fontsize',8)

f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_Antibody.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_Antibody.pdf']),'ContentType','vector',...
            'BackgroundColor','none');