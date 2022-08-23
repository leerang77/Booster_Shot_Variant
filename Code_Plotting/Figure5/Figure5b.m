%%
load(fullfile('..','Results',[mfilename, '_Epitope_Masking_Steric_Hindrance_Agconc.mat']))
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
f = figure;
t = tiledlayout(1,1,'Padding','tight');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 5.2, 4.7];
nexttile;
Linestyles = {'-','--',':','-.'};

for ep=1:2
h(ep)=plot(param{2}.tspan_summary, Ceff(ep,:), 'color', colors{ep}, 'Linewidth', 1, 'Linestyle', Linestyles{ep});
hold on
end

set(gca,'YScale','log')
xlim([0,30])
xticks(0:15:30)
xlabel('Time after Vax2 (day)', 'fontsize', 8)
ylabel('Effective Ag concentration (nM)', 'fontsize', 8)
ylim(10.^([-4,1]))
yticks(logspace(-5, 1, 4))
leg = legend({'Dominant', 'Subdominant'}, 'fontsize', 8, 'Location', 'Northeast');
leg.ItemTokenSize = [15,5];
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_Epitope_Masking_Steric_Hindrance_Agconc.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_Epitope_Masking_Steric_Hindrance_Agconc.pdf']),'ContentType','vector',...
            'BackgroundColor','none');