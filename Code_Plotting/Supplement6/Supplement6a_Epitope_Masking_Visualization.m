fprintf(mfilename)
Kds = [10^-8, 10^-9, 10^-10];
x = logspace(-9.5, -6, 20); % M
x_mass = x*150000*10^6/1000; %ug/ml
cmap=[44,162,95; 65,182,196; 129,15,124]/256;
figure
t = tiledlayout(1,1,'Padding','compact');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 5, 6.8];
nexttile;
p=[];
for Kd=Kds
    p(Kd==Kds)=loglog(x_mass,Kd./(Kd+x),'color',cmap(Kd==Kds,:), 'LineWidth', 1);
    names{Kd==Kds} = sprintf('K_d = 10^{%d} M', -log10(Kd));
    hold on
end
vax = zeros(2,2);
vax(1,:) = [3.6/5, 1/(1.06*10^9)]; %amount in ug/ml, and Kd
vax(2,:) = [25.5/5, 1/(3.6*10^9)];
free = vax(:,2)./(vax(:,2)+vax(:,1)*1000/(150000*10^6));
p(4)=scatter(vax(1,1), free(1), 50, [164,164,164]/256, 'd',  'filled');
p(5)=scatter(vax(2,1), free(2), 50, [10,10,10]/256, 'd',  'filled');

set(gca,'FontSize',8)
leg = legend(p, [names,{'~3 weeks Vax1','~2 weeks Vax2'}],'location','northoutside','numcolumns',2, 'fontsize', 6);
leg.ItemTokenSize = [8,5];
xticks([0.1, 1, 10, 100])
xlabel(t, 'Antibody concentration (\mug/ml)', 'fontsize', 8)
ylabel(t, 'Fraction of epitope remaining unmasked', 'fontsize', 8)
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'.pdf']),'ContentType','vector',...
            'BackgroundColor','none');

% figure
% for Kd=Kds
%     loglog(x,Kd./(Kd+x),'DisplayName', ['Kd = ', num2str(Kd), ' M'])
%     hold on
% end
% set(gca,'FontSize',16)
% xlabel('Ab conc (M)')
% legend('show','location','best')
% ylabel('Fraction of free epitope')
% set(gcf,'Position',[100,100,400,350]);

% figure
% x = logspace(-3, 2, 20);
% loglog(x, 1./(x+1))
% set(gca,'FontSize', 14);
% xlabel('[Ab]*Ka')
% ylabel('Fraction of free epitope')
% set(gcf,'Position',[100,100,400,350]);