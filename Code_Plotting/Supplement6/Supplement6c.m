%%
addpath('..');
load(fullfile('..','Results',[mfilename, '_Epitope_Masking_Steric_Hindrance_GC.mat']));
rng(20);
nn = randsample(M_GC, 100)';
t = tiledlayout(1,1,'Padding','compact');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 5.3, 4.8];
for i=2:2
    nexttile;
    p = plotGCNum(param{i}, result{i}, 2, [], colors);
end
h = findobj(gca,'Type','line');
% h(1).LineWidth = 1;
% h(2).LineWidth = 1;
ylim([1,10^4.3])
yticks(logspace(0,4,5))
leg = legend(p(1:2), 'Dominant', 'Subdominant');
leg.ItemTokenSize = [10,5];
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_Epitope_Masking_Steric_Hindrance_GC.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_Epitope_Masking_Steric_Hindrance_GC.pdf']),'ContentType','vector',...
            'BackgroundColor','none');