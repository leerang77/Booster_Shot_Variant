%% Vax1
load(fullfile('..','Results',[mfilename, '_Vax1Scan_E1h_dE12.mat']));
figure
subplot(1,2,1)
h1 = heatmap(x,y,log10(memnum{2}'),'CellLabelColor','none', 'Colormap', flipud(hot));
h1.ColorbarVisible='off';
xlabel('E_1^h')     
ylabel('dE_{12}')
set(gcf,'Position',[100,100,360,170])
caxis([2,5]);
c.TickLabels = {'10^2','10^3','10^4','10^5'};
title('Vax1 1m GC')
set(gca,'Fontsize',8)

%% Vax2
load(fullfile('..','Results',[mfilename, '_Vax2Scan_E1h_dE12.mat']));
subplot(1,2,2)
heatmap(x,y,log10(memnum{2}'),'CellLabelColor','none', 'Colormap', flipud(hot))
xlabel('E_1^h')     
ylabel('dE_{12}')
caxis([2,5]);
title('Vax2 5m GC')
set(gca,'Fontsize',8)
sgtitle('Memory cells', 'fontsize', 10, 'fontweight', 'bold')
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_Scan_E1h_dE12.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_Scan_E1h_dE12.pdf']),'ContentType','vector',...
            'BackgroundColor','none');