fprintf(mfilename)
addpath(fullfile('..','Code_Parameter_Generation'));
addpath(fullfile('..','Code_Result_Analysis'));
% Load prime result, no Ep masking
p = base_case_parameters();
p.vaxnum = 2; p.tmax = 180;
result = cell(1,3); param = cell(1,3);


pc = cell(3,2);
mem = cell(3,2);
outputprob = [0.03, 0.05, 0.10];
outputpcfrac = [0.1, 0.5];
figure;
t = tiledlayout(1,1,'Padding','compact');
colors = cell(3,1);
colors{1} = {'#feb24c', '#feb24c'};
colors{2} = {'#7fcdbb', '#2c7fb8'};
colors{3} = {'#c994c7', '#dd1c77'};

% colors = {[0 0.4470 0.7410],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980]};
for i=1:3
    for j=1:2
        pathPrefix = fullfile('outputChange',sprintf('outputprob_%.2f',outputprob(i)),sprintf('outputpcfrac_%.2f', outputpcfrac(j)));
        result{i,j} = load_result(vaxnum, T, k, numshot, E1h, dE12, p, n_ep, ...
            masking, production, exfpb, k_deposit, IgM0, Ag0, C0, Ageff,...
            w1, w2, tmax, first, last, steric, memToGC, pathPrefix);
        param{i,j} = result{i,j}.param;
        [pc{i,j},mem{i,j}] = MemPCnum(param{i,j}, result{i,j});
        semilogy(param{i,j}.tspan_summary, mem{i,j}{1}(1,:), '--', 'color', colors{i}{j});
        hold on
        semilogy(param{i,j}.tspan_summary, mem{i,j}{1}(2,:), 'color', colors{i}{j});
    end
    legend(arrayfun(@(x) num2str(x), outputpcfrac, 'UniformOutput', false))
    title(num2str(outputprob(i)));
end
save('Supplement4E.mat','param','pc', 'mem', 'outputprob', 'outputpcfrac', 'colors')
f = gcf;
exportgraphics(f,fullfile('figures','Supplement_ChangingOutput_MemPC.pdf'),'ContentType','vector',...
            'BackgroundColor','none')

%%
load('Supplement4E.mat')
outputprob = [0.03, 0.05, 0.10];
outputpcfrac = [0.1, 0.5];
colors{1} = {'#feb24c', '#feb24c'};
colors{2} = {'#7fcdbb', '#2c7fb8'};
colors{3} = {'#c994c7', '#dd1c77'};

figure
t = tiledlayout(1,1,'Padding','compact');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 7, 5];
labels = cell(1,length(outputprob)*length(outputpcfrac));
idx = 1;
linestyle = {'-','--'};
nexttile;
semilogy(param{i,j}.tspan_summary, mem{i,j}{1}(1,:), '-', 'color', 'k');
hold on
for i=1:3
    for j=1:2
        semilogy(param{i,j}.tspan_summary, mem{i,j}{1}(2,:), 'color', colors{i}{2}, 'LineStyle', linestyle{j});
        labels{idx} = sprintf('p_1=%.2f, p_2=%.1f', outputprob(i), outputpcfrac(j));
        idx = idx+1;
    end
end
ylim([10^2, 10^6.5])
xlim([0,150])
legend([{'EGC'},labels], 'location', 'southeast', 'fontsize', 6)

xlabel('Time (day)', 'fontsize', 8)
ylabel('Antibody Titer', 'fontsize', 8)
% leg.ItemTokenSize = [15,5];
    f = gcf;
exportgraphics(f,fullfile('figures','Supplement4E.pdf'),'ContentType','vector',...
        'BackgroundColor','none');