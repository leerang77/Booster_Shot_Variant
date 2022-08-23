suffixes = {'','_EpMasking'};
for i=1:2
    suffix = suffixes{i};
%%
load(fullfile('..','Results',[mfilename, '_Comparisons', suffix,'.mat']))
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
figure
t = tiledlayout(1,1,'Padding','compact');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 8.2, 5];
nexttile;
x = 1:4;
scatter(x,WT_titer(1,[1,2,4,3]),20,colors{1},'filled')
hold on
scatter(x,WT_titer(2,[1,2,4,3]),20,colors{2},'filled')
scatter(x,Var_titer(1,[1,2,4,3]),20,colors{1})
scatter(x,Var_titer(2,[1,2,4,3]),20,colors{2})
set(gca,'YScale','log')
box on
leg = legend('WT-Dominant', 'WT-Subdominant', 'Variant-Dominant', 'Variant-Subdominant', 'location','northoutside','fontsize', 6);
leg.ItemTokenSize = [15,5];
leg.Layout.Tile = 'east';
xlim([0,5])
xticks([1,2,3,4,5])
xticklabels({'Vax2','Vax3','Vax4','Vax3-Short'})
ylabel('Titer (C_{Ab}K_{a})','fontsize',8)
f = gcf;
savefig(f, fullfile('..','figures',['Figure6_Comparisons',suffix, '.fig']))
exportgraphics(f,fullfile('..','figures',['Figure6_Comparisons',suffix, '.pdf']),'ContentType','vector',...
            'BackgroundColor','none')
%%
figure
t = tiledlayout(1,1,'Padding','compact');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 7, 6];
nexttile;
x=1:4;
scatter(x,sum(WT_titer(:,[1,2,4,3])),20,'MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0])
hold on
scatter(x,sum(Var_titer(:,[1,2,4,3])),20,'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',[0,0,0])
set(gca,'YScale','log')
leg = legend('WT', 'Variant', 'location','southeast','fontsize', 8);
leg.ItemTokenSize = [15,5];
% leg.Layout.Tile = 'east';
xlim([0,5])
ylim([10^-2,10^4])
xticks([1,2,3,4,5])
xticklabels({'Vax2 1.3m','Vax3 1m','Vax4 1m','Vax3-Short 1m'})
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',8)
ylabel('Titer (C_{Ab}K_{a})','fontsize',8)
box on
f = gcf;
savefig(f, fullfile('..','figures',['Figure6_Comparisons_sum',suffix, '.fig']))
exportgraphics(f,fullfile('..','figures',['Figure6_Comparisons_sum',suffix, '.pdf']),'ContentType','vector',...
            'BackgroundColor','white')
end