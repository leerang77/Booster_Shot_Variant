%%
load(fullfile('..','Results',[mfilename, '_PCHistogram.mat']));
plotPCHistogram(ctrs, affcnts_gcpc, affcnts_egcpc)
%%
function plotPCHistogram(ctrs, affcnts_gcpc, affcnts_egcpc)
num_sim=10;
barcolors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.3010 0.7450 0.9330], [0.9290 0.6940 0.1250]};
strainName = {'WT', 'Variant'};
for ep=1:2
    figure(ep)
    t(ep) = tiledlayout(1,2,'Padding','tight');
    t(ep).Units = 'centimeters';
    t(ep).OuterPosition = [3, 3, 4.5, 5];
end
ax = gobjects(2,2);
for strain=1:2    
    figure(1);
    ax(1,strain) = nexttile;
    b = bar(ctrs, [ affcnts_gcpc{2}{strain,1}; affcnts_gcpc{3}{strain,1};affcnts_egcpc{2}{strain,1};]/num_sim, 'stacked', 'FaceColor', 'flat', 'EdgeColor', 'flat');
    b(1).CData = [195,223,247]/256;
    b(2).CData = barcolors{1};
    b(3).CData = barcolors{3};
    title({strainName{strain},''})
    figure(2);
    ax(2,strain) = nexttile;
    b = bar(ctrs, [affcnts_gcpc{2}{strain,2}; affcnts_gcpc{3}{strain,2};affcnts_egcpc{2}{strain,2};]/num_sim, 'stacked', 'FaceColor', 'flat', 'EdgeColor', 'flat');
    b(1).CData = [247,206,195]/256;
    b(2).CData = barcolors{2};
    b(3).CData = barcolors{4};    
    title({strainName{strain},''})
end
names = {'dominant', 'subdominant'};
    for j=1:2
ylim(ax(1,j),[0,3*10^5])
ylim(ax(2,j),[0,1.5*10^4])
    end
ylim(ax(1,2),[0,1.5*10^4])   
ax(2,1).YAxis.Exponent = 4;
ax(1,2).YAxis.Exponent = 4;
ax(2,2).YAxis.Exponent = 4;
for ep=1:2
    figure(ep)
    leg = legend({'1.3m GC', '5m GC', 'EGC'}, 'fontsize', 6, 'Orientation','Horizontal', 'NumColumns', 3, 'Location', 'NorthOutside');
    leg.ItemTokenSize = [10,3];
    leg.Layout.Tile = 'North';
    xlabel(t(ep), ['Binding affinity (-log_{10}K_d)'], 'fontsize', 8)
    ylabel(t(ep), ['Vax2 ', names{ep},' plasma cells'], 'fontsize', 8)
    f = gcf;
    savefig(f, fullfile('..','figures',['Figure2f_PCHist_',strainName{strain},'.fig']))
    exportgraphics(f,fullfile('..','figures',['Figure2f_PCHist_',names{ep},'.pdf']),'ContentType','vector',...
            'BackgroundColor','none');
end


end
