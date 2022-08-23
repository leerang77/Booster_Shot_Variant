%%
load(fullfile('..','Results',[mfilename, '_NoMaskingvsMasking.mat']))
f = figure;
set(gcf,'Position',[100,100,1000,500])
t = tiledlayout(1,1,'Padding','compact');
t.Units = 'centimeters';
t.OuterPosition = [2, 2, 5.4, 4.7];
a1 = nexttile;
x_indx = [1,2];
p = bar(x_indx, [sum(vax3.WT_titer,1);sum(vax3.Var_titer,1)]','grouped', 'FaceColor', 'flat');
colors = [0, 0, 0;
          1,1,1];
p(1).FaceColor = colors(1,:);
p(2).FaceColor = colors(2,:);
set(gca,'YScale','log')
leg = legend(p(1:2), {'WT', 'Variant'},'fontsize', 6, 'location', 'north', 'orientation', 'horizontal');
leg.ItemTokenSize = [15,5];
% leg.Layout.Tile = 'East';
xlim([0.25,2.75])
ylim(10.^([0,5]))
yticks(logspace(-1,5,4))
ylabel('Titer (C_{Ab}K_{a})','fontsize',8)
xticks([1,2])
xticklabels({'No Masking','Masking'})
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_NoMaskingvsMasking.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_NoMaskingvsMasking.pdf']),'ContentType','vector',...
            'BackgroundColor','none');