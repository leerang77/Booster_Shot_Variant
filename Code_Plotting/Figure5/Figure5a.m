fprintf(mfilename)
data = [0	0	0.023
0	0	0.037
0.04	0.117	0.184
0.179	0.179	0.155
0	0.075	0.177
0.258	0.147	0.103
0.135	0.127	0.0809
0	0.02	0
0.04	0	0
0.061	0.059	0.0376
0.287	0.27	0.2];

blocking = zeros(4,11);
mapping = {1, 2, 3, 4, [1,4], [1,2], [2,3], [3,4], [1,2,3], [1,2,4], []};
for i=1:4
   for j=1:10
      if ismember(i, mapping{j})
          blocking(i,j) = 1;
      end
   end
end

sites_interfering = blocking*data;
c = [236,226,240;
    166,189,219;
    28,144,153]/256;
figure
t = tiledlayout(1,1,'Padding','compact');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 5.1, 4.9];
nexttile;
h = bar([1,2,3,4],sites_interfering,'FaceColor','flat');
for i=1:3
    h(i).CData = c(i,:);
end
ylim([0,0.55])
xticklabels({'1','2','3','1/4'})
xlabel('Reference Antibody Class', 'fontsize', 8)
ylabel('Fraction of serum antibodies blocked', 'fontsize', 8)
leg = legend('Vax2 1.3m', 'Vax2 5m', 'Vax3 1m');
leg.ItemTokenSize = [8,5];
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'.pdf']),'ContentType','vector',...
            'BackgroundColor','none');