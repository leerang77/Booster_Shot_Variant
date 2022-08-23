%%
load(fullfile('..','Results',[mfilename, '_Unique_Sequences_0.mat']))

% Plot
figure
t = tiledlayout(1,1,'Padding','tight');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 7, 7];
ax = nexttile;
colors = {[117,107,177]/256,[253,174,97]/256,[28,144,153]/256};

plot(1:length(GC_seq_sizes{1}), GC_seq_sizes{1}, 'Linewidth', 1.5, 'Color', colors{1})
hold on
plot(1:length(EGC_seq_sizes{2}), EGC_seq_sizes{2}, 'Linewidth', 1.5, 'Color', colors{3});
plot(1:length(GC_seq_sizes{2}), GC_seq_sizes{2}, 'Linewidth', 1.5, 'Color', colors{2});

colormap(ax,'turbo')
% set(ax, 'YScale', 'log')
leg = legend({'Vax1 1m GC', 'Vax2 EGC', 'Vax2 5m GC'}, 'Location', 'NorthEast');
leg.ItemTokenSize = [15, 5];
xlabel('Rank of sequence')
ylabel('Number of B cells with identical sequence')
set(gca,'YScale','log')
xlim([0,300])

f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_Unique_Sequences.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_Unique_Sequences.pdf']),'ContentType','vector',...
            'BackgroundColor','none');