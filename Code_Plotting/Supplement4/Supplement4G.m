%%
load(fullfile('..','Results',[mfilename, '_ChangingOutput.mat']))

%Plot Vax1-2 antibodies
linestyles = {'-','--','--'};
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
memColors = {[253,174,97]/256, [28,144,153]/256};

figure
set(gcf,'Position',[100,100,1000,500])
t = tiledlayout(1,3,'Padding','tight');
t.Units = 'centimeters';
t.InnerPosition = [3, 3, 11, 2.5];
p = [];
nexttile;
for i=[1,3]
    idx = (1+4*28):4:(1+178*4);
    p(i) = plot(time(idx)-28,WT_titer{i}(idx), linestyles{i}, 'color', colors{1}, 'LineWidth', 1)
    hold on
    p(i+1) = plot(time(idx)-28,Variant_titer{i}(idx), linestyles{i}, 'color', colors{2}, 'LineWidth', 1)
end
set(gca,'YScale','log')
ylim([10^-1,10^2.5])
yticks(logspace(-1,3,3))

set(gca,'fontname','arial')
xlabel('Time after Vax2 (day)','fontsize',8)
ylabel('Titer','fontsize',8)
leg = legend(p([1,3,2,4]), {'0.1','0.5', '0.1', '0.5'},  'numcolumns', 2, 'location', 'northoutside', 'fontsize', 6);
leg.ItemTokenSize = [18,3];
title('Changing p_2')
%Plot Vax1-2 Memory number
pl = [];
nexttile;
for i=[1,2]
    idx = (1+4*28):4:(1+178*4);
    pl(i*2-1) = plot(time(idx)-28,sum([memnum{i,1}(1,idx); memnum{i,2}(1,idx)]), linestyles{i}, 'color', memColors{1}, 'LineWidth', 1)
    hold on
    pl(i*2) = plot(time(idx)-28,sum([memnum{i,1}(2,idx); memnum{i,2}(2,idx)]), linestyles{i}, 'color', memColors{2}, 'LineWidth', 1)
end
ylim([10^5,10^6])
set(gca,'fontname','arial')
xlabel('Time after Vax2 (day)','fontsize',8)
ylabel('Memory cells','fontsize',8)
set(gca,'YScale','log')
leg = legend(pl([1,3,2,4]), {'0.05','0.1', '0.05', '0.1'}, 'numcolumns', 2, 'location', 'northoutside', 'fontsize', 6);
leg.ItemTokenSize = [18,3];
title('Changing p_1')
%Plot Antibody Comparison
nexttile;
x_indx = [1,3];
pl = [];
markerstyles = {'o','*','x'};
for i=1:3
% pl(1) = plot(x_indx,titers{1}(i,:), [linestyles{i},'o'], 'MarkerSize', 3, 'color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
% hold on
% pl(2) = plot(x_indx,titers{2}(i,:), [linestyles{i},'o'], 'MarkerSize', 3, 'color', 'k', 'MarkerFaceColor', [1,1,1], 'MarkerEdgeColor', 'k');
pl(i*2-1) = plot(x_indx,titers{1}(i,[1,3]), ['-',markerstyles{i}], 'LineWidth', 0.2, 'MarkerSize', 7, 'color', 'k', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', colors{1});
hold on
pl(i*2) = plot(x_indx,titers{2}(i,[1,3]), ['-',markerstyles{i}], 'LineWidth', 0.2, 'MarkerSize', 7, 'color', 'k', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', colors{2});
end
set(gca,'YScale','log')
leg = legend(pl([1,3,5,2,4,6]),{'0.05, 0.1','0.05, 0.5','0.1, 0.1','0.05, 0.1','0.05, 0.5','0.1, 0.1'}, 'numcolumns', 2, 'location', 'northoutside', 'fontsize', 6);
leg.ItemTokenSize = [18,3];
set(gca,'YScale','log')
ylabel('Titer','fontsize',8)
xticks([1,3])
xlim([0,4])
ylim([10^-1.5, 10^4])
yticks(logspace(-1,3,3))
% yticklabels(cellfun(@(x) num2str(x), num2cell(logspace(-1,4,6)), 'UniformOutput', false))
xticklabels({'Vax2 1.3m', 'Vax3 1m'})
title('Changing p_1, p_2')
set(gca, 'xticklabel', get(gca, 'xticklabel'), 'fontsize', 8)
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_ChangingOutput.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_ChangingOutput.pdf']),'ContentType','vector',...
            'BackgroundColor','none');