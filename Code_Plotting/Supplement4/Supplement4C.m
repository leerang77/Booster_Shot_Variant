%% Vax1
load(fullfile('..','Results',[mfilename, '_Vax1Scan_dE12_p.mat']));
figure
subplot(1,2,1)
h1 = heatmap(x,y,log10(memnum{2}'),'CellLabelColor','none', 'Colormap', flipud(hot));
h1.ColorbarVisible = 'off';
xlabel('dE_{12}')
ylabel('p')
caxis([2,5]);
title('Vax1 1m GC')

%% Vax2
load(fullfile('..','Results',[mfilename, '_Vax2Scan_dE12_p.mat']));
subplot(1,2,2)
heatmap(x,y,log10(memnum{2}'),'CellLabelColor','none', 'Colormap', flipud(hot))
xlabel('dE_{12}')
ylabel('p')
caxis([2,5]);
title('Vax2 5m GC')
set(gcf,'Position',[100,100,360,170])
set(gca,'Fontsize',8)
sgtitle('Memory cells', 'fontsize', 10, 'fontweight', 'bold')
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_Scan_dE12_p.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_Scan_dE12_p.pdf']),'ContentType','vector',...
            'BackgroundColor','none');