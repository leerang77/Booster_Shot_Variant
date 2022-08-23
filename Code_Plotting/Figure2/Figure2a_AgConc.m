fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
%%
[result, param] = loadFigure2Data();
result{2} = result{3}; param{2} = param{3};
agconc_mean = cell(1,2);
for i=1:2
   agconc_mean{i} = PlotAgConc(param{i}, result{i}); 
end
agconc_x = [param{1}.tspan_summary, param{2}.tspan_summary(2:end)+28];
agconc_y = [agconc_mean{1}, agconc_mean{2}(:,2:end)];
clear result
outputFile = fullfile('..','Results',[mfilename, '.mat']);
save(outputFile)

%% Ag concentration
function agconc_mean = PlotAgConc(param, result)
agconc_mean = zeros(2,length(param.tspan_summary));
figure
for i=1:length(result.conc)
    agconc = squeeze(result.conc{i}.concarray(1,:,:));
    agconc = [agconc(1,:); sum(agconc(2:end, :),1)];
    
    agconc_mean = agconc_mean + agconc;
    agconc_mean(isnan(agconc_mean)) = 0;
    if i==1
    p1 = semilogy(param.tspan_summary, agconc(1,:), 'color', [0, 0, 0, 1]);
    hold on
    p2 = semilogy(param.tspan_summary, sum(agconc(2,:),1), 'color', [0, 0, 0, 1], 'LineStyle', '--');
    end
end
agconc_mean = agconc_mean/length(result.conc);
legend([p1, p2], {'soluble', 'IC-FDC'})
xlabel('Time(months)')
ylabel('Ag concentration (nM)')
end
