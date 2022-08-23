function [result, param] = loadFigure2Data()
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultTextFontSize', 8);
p = base_case_parameters();
result = cell(1,3); param = cell(1,3);
[result{1},dirname] = load_result(p);
param{1} = result{1}.param;

p.vaxnum = 2;
p.tmax = 42;
[result{2},dirname] = load_result(p);
param{2} = result{2}.param;

p.vaxnum = 2;
p.tmax = 180;
[result{3},dirname] = load_result(p);
param{3} = result{3}.param;

M_GC = 2000;
colors = {[0 0.4470 0.7410],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980]};
colors = {colors{1}, colors{2:param{1}.n_ep-1}, colors{end}};
end