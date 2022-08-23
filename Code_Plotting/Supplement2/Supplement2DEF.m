%%
load(fullfile('..','Results',[mfilename, '_MemoryScatter.mat']))
plotMemoryScatter(memaff, memaff_egc)

%%
function plotMemoryScatter(memaff, memaff_egc)
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
h = cell(1,3);
h{1} = plotScatter(memaff{1}, colors);
title(h{1},'Vax1 1m GC-derived memory cells','fontsize',10)
f = gcf;
exportgraphics(f,fullfile('..','figures','Supplement2_Memscatter1.png'),'Resolution', 600,...
    'BackgroundColor','none');
h{2} = plotScatter(memaff{3}, colors);
title(h{2}, 'Vax2 5m GC-derived memory cells','fontsize',10)
f = gcf;
exportgraphics(f,fullfile('..','figures','Supplement2_Memscatter2.png'),'Resolution', 600,...
    'BackgroundColor','none');
h{3} = plotScatter(memaff_egc{3}, colors);
title(h{3}, 'Vax2 EGC-derived memory cells','fontsize',10)
f = gcf;
exportgraphics(f,fullfile('..','figures','Supplement2_Memscatter3.png'),'Resolution', 600,...
    'BackgroundColor','none');
end

function t = plotScatter(affs, colors)
figure;
t = tiledlayout(1,2,'Padding','Tight');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 8, 4.3];
nexttile;
scatterUnique(affs{1,1},affs{2,1},colors{1});
xlabel('WT affinity')
ylabel('Variant affinity')
ylim([5, 12])
xlim([5, 12])
xticks(6:2:12)
xticklabels([{'<6'},cellfun(@(x) num2str(x), num2cell(8:2:12), 'UniformOutput', false)])
yticks(6:2:12)
yticklabels([{'<6'},cellfun(@(x) num2str(x), num2cell(8:2:12), 'UniformOutput', false)]) 
nexttile;
scatterUnique(affs{1,2},affs{2,2},colors{2});
xlabel('WT affinity')
ylabel('Variant affinity')
ylim([5, 12])
xlim([5, 12])
xticks(6:2:12)
xticklabels([{'<6'},cellfun(@(x) num2str(x), num2cell(8:2:12), 'UniformOutput', false)])
yticks(6:2:12)
yticklabels([{'<6'},cellfun(@(x) num2str(x), num2cell(8:2:12), 'UniformOutput', false)])
end

function scatterUnique(x, y, color)
 A = [x;y];
 [uniqueA,~,Idx] = unique(A','rows');
 cnts = histcounts(Idx, 1:max(Idx)+1);
 size = min(0.004*cnts,50);
 uniqueA = uniqueA';
 scatter(uniqueA(1,:), uniqueA(2,:), size, color, 'filled')
end