%%
addpath('..');
load(fullfile('..','Results',[mfilename, '_Epitope_Masking_Steric_Hindrance_GC.mat']));
rng(20);
nn = randsample(M_GC, 100)';
t = tiledlayout(1,1,'Padding','compact');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 5.5, 5];
for i=2:2
    nexttile;
    p = plotGCNum(param{i}, result{i}, 2, [], colors);
end
xticks([0,60,120])
ylim([1,10^4])
yticks(logspace(0,4,3))
leg = legend(p(1:2), 'Dominant', 'Subdominant');
leg.ItemTokenSize = [10,5];
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_Epitope_Masking_Steric_Hindrance_GC.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_Epitope_Masking_Steric_Hindrance_GC.pdf']),'ContentType','vector',...
            'BackgroundColor','none');