function p = plotGCNum(param, result, vaxnum, nn, colors)
xlength = length(param.tspan_summary);
for i=nn
    for ep=1:param.n_ep
    p(ep)=semilogy(param.tspan_summary, result.gc.numbytime(i,xlength*(ep-1)+1:xlength*ep,1),'-','Color',[colors{ep}, 0.2]);
    hold on
    end
end

% legend(p(:),legendCellEp)
totalnum = squeeze(sum(result.gc.numbytime, 1))/size(result.gc.numbytime,1);
% yyaxis right
for ep=1:param.n_ep
    p(ep)=semilogy(param.tspan_summary, totalnum(xlength*(ep-1)+1:xlength*ep,1),'-', 'LineWidth', 3, 'Color',[colors{ep}]);
    hold on
end
set(gca,'fontname','arial')
xlabel(sprintf('Time after Vax%d (day)', vaxnum),'fontsize',8)
ylabel('GC B cells','fontsize',8)
ylim([1,5000])
xlim([0, max(param.tspan_summary)])
if max(param.tspan_summary)==28
    xticks(0:15:30)
    xlim([0,30])
elseif max(param.tspan_summary)==180
    xticks(0:30:150)
    xlim([0,150])
end
end