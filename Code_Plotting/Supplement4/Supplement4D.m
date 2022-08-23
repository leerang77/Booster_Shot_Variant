%%
load(fullfile('..','Results',[mfilename, '_Changingw2.mat']))
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
figure
t = tiledlayout(1,2,'Padding','tight');
t.Units = 'centimeters';
t.InnerPosition = [3, 3, 6, 2.5];

p = gobjects(1,4);
for i=1:2
    memnum{i}(memnum{i}<10) = 10;
end
nexttile;
p(1) = plot(w2s,memnum{1}(1,:), '-o','MarkerSize', 4, 'color', colors{1});
hold on
p(3) = plot(w2s,memnum{1}(2,:), '-o','MarkerSize', 4, 'color', colors{2});
set(gca,'YScale','log')
ylim([10^1, 10^5])
yticks(logspace(1,5,3))
xticks([.3, 1])
xlabel('K', 'fontsize',8)
title('Vax1 1m GC')

nexttile;
p(2) = plot(w2s,memnum{2}(1,:), '-o','MarkerSize', 4, 'color', colors{1});
hold on
p(4) = plot(w2s,memnum{2}(2,:), '-o','MarkerSize', 4, 'color', colors{2});
ylim([10^2, 10^6])
xticks([.3, 1])

yticks(logspace(2,6,3))
set(gca,'Yscale','log')
title('Vax2 5m GC')

ylabel(t,'Memory cells','fontsize',8)
xlabel('K', 'fontsize',8)
leg = legend(p([2,4]),{'Dominant', 'Subdominant'}, 'Location', 'Southeast', 'fontsize', 6);
leg.ItemTokenSize = [10,3];
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_Changingw2.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_Changingw2.pdf']),'ContentType','vector',...
            'BackgroundColor','none');
