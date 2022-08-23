%% GC Vax1
% load(fullfile('..','Data','Data_Prime','1_0_0_1_7_0.4_0.1_2_0_1_1_24_0.01_10_0.008_0.01_1_0.5_28','1_to_200.mat'));
load('..\Data_Prime\1_0_0_1_7_0.4_0.1_2_0_0.2_0_24_0.01_10_0.008_0.01_1_0.5_28\1_to_200.mat');
f = 7;
output = result.output.finalmem;
naiveaff = squeeze(result.naive(:,:,3));
n_GC = size(output,1);
num_uniques = [];
for k=1:n_GC
   mem = output(k,output(k,:,3)>f,:);
   if ~isempty(mem)
   uniqueProperty = (mem(1,:,1).^2+mem(1,:,3)+mem(1,:,4)+mem(1,:,5));
   nums = histcounts(uniqueProperty, [unique(uniqueProperty),inf]);
   num_uniques = [num_uniques, nums];
   end
end
GC{1} = sort(num_uniques,'descend');
%% GC Vax2
% load(fullfile('..','Data','Data_Secondary','2_0_0_1_7_0.4_0.1_2_0_1_1_24_0.01_10_0.008_0.01_1_0.5_180','1_to_200.mat'));
load('..\Data_Secondary\2_0_0_1_7_0.4_0.1_2_0_0.2_0_24_0.01_10_0.008_0.01_1_0.5_180\1_to_200.mat');
% load('..\Data_Secondary\2_0_0_1_7_0.4_0.1_2_0_1_1_24_0.01_10_0.008_0.01_1_1_6\1_to_200.mat');
f = 7;
output = result.output.finalmem;
naiveaff = squeeze(result.naive(:,:,3));
n_GC = size(output,1);
num_uniques = [];
for k=1:n_GC
   mem = output(k,output(k,:,3)>f,:);
   if ~isempty(mem)
   uniqueProperty = (mem(1,:,1).^2+mem(1,:,3)+mem(1,:,4)+mem(1,:,5));
   nums = histcounts(uniqueProperty, [unique(uniqueProperty),inf]);
   num_uniques = [num_uniques, nums];
   end
end
GC{2} = sort(num_uniques,'descend');
%% EGC
uniqueidx = result.memoryCellsEGC(9,:)+result.memoryCellsEGC(4,:);
uniqueidx = uniqueidx(result.memoryCellsEGC(3,:)>f);
cnts = arrayfun(@(x) sum(x==uniqueidx), unique(uniqueidx));
EGC{2} = sort(cnts, 'descend');

%% Plot
figure
set(gca,'YScale','log')
xlabel('rank')
ylabel('size')
set(gcf,'Position',[100,100,300,300])
figure
t = tiledlayout(1,1,'Padding','tight');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 5, 5];
ax = nexttile;
colors = {'#ebac23','#006e00','#00b25d'};
plot(1:length(GC{1}), GC{1}, 'Linewidth', 1.5, 'Color', colors{1})
hold on
plot(1:length(EGC{2}), EGC{2}, 'Linewidth', 1.5, 'Color', colors{3});
plot(1:length(GC{2}), GC{2}, 'Linewidth', 1.5, 'Color', colors{2});

colormap(ax,'turbo')
set(ax, 'YScale', 'log')
leg = legend({'Vax1 1m GC', 'Vax2 5m EGC', 'Vax2 5m GC'}, 'Location', 'NorthEast');
leg.ItemTokenSize = [15, 5];
xlabel('rank')
ylabel('size')
xlim([0,1000])
yticks(logspace(0,5,6));
f = gcf;
exportgraphics(f,fullfile('figures','Supplement_Unique_Sequences_Size.pdf'),'ContentType','vector',...
    'BackgroundColor','none');
