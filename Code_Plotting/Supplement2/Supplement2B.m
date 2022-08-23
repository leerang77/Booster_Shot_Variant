%%
load(fullfile('..','Results',[mfilename, '_14dGCDiversity.mat']))
colors = {[0 0.4470 0.7410],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980]};
colors = {colors{1}, colors{2:param{1}.n_ep-1}, colors{end}};
figure
t = tiledlayout(1,1,'Padding','tight');
t.Units = 'centimeters';
t.InnerPosition = [3, 3, 3.7, 3.7];
nexttile;
b = bar(ctrs, clonality_hist, 1, 'FaceColor', 'flat', 'EdgeColor', 'flat');
b(1).CData = [117,107,177]/256;
b(2).CData = [28,144,153]/256;
legend({'Vax1','Vax2'})
xlabel('Occupancy of the largest lineage')
ylabel('Fraction of simulated GCs')
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_14dGCDiversity.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_14dGCDiversity.pdf']),'ContentType','vector',...
            'BackgroundColor','none');