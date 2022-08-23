%%
load(fullfile('..','Results',[mfilename, '_Naive_Entry.mat']))
colors = {[0 0.4470 0.7410],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980]};
colors = {colors{1}, colors{2:param{1}.n_ep-1}, colors{end}};
figure
t = tiledlayout(1,1,'Padding','tight');
t.Units = 'centimeters';
t.InnerPosition = [3, 3, 3.7, 3.7];
nexttile;
tspan = 0:length(naive_entry{1}(1,:))-1;
plot(tspan, naive_entry{1}(1,:), '-o', 'Color', colors{1})
hold on
plot(tspan, naive_entry{2}(1,:), '-o','Color', colors{2})
plot(tspan, naive_entry{1}(2,:), '--x', 'Color', colors{1})
plot(tspan, naive_entry{2}(2,:), '--x', 'Color', colors{2})
xlim([0,8])
xlabel('Time (day)')
ylabel({'Number of naive B cells' 'that entered GC'})
% xlim([0,8])
ylim([0,110])

leg = legend({'Vax1 Dominant', 'Vax1 Subdominant', 'Vax2 Dominant', 'Vax2 Subdominant'}, 'fontsize', 6, 'location', 'northwest');
leg.ItemTokenSize = [15,5];
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_Naive_Entry.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_Naive_Entry.pdf']),'ContentType','vector',...
            'BackgroundColor','none');


figure
t = tiledlayout(1,1,'Padding','tight');
t.Units = 'centimeters';
t.InnerPosition = [3, 3, 3.7, 3.7];
nexttile;
plot(affs, squeeze(num_by_aff(1,1,:)), '-o', 'Color', colors{1})
hold on
plot(affs, squeeze(num_by_aff(2,1,:)), '-o', 'Color', colors{2})
plot(affs, squeeze(num_by_aff(1,2,:)), '-x', 'Color', colors{1})
plot(affs, squeeze(num_by_aff(2,2,:)), '-x', 'Color', colors{2})
set(gca,'Yscale','log')

leg = legend({'Vax1 Dominant', 'Vax1 Subdominant', 'Vax2 Dominant', 'Vax2 Subdominant'}, 'fontsize', 6, 'location', 'northeast');
leg.ItemTokenSize = [15,5];
xlabel('Germline affinity (-log_{10}K_d)')
ylabel({'Number of naive B cells' 'that entered GC'})
ylim([10^-4, 10^4])
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'_Naive_Entry.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'_Num_By_Aff.pdf']),'ContentType','vector',...
            'BackgroundColor','none');
