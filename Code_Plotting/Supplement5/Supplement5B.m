load(fullfile('..','Results','Supplement5A_Unique_Sequences_0.mat'))

%% Calculate 
Sarray = 10:5:100;
EGCTPRs = zeros(size(Sarray));
EGCTPR2s = zeros(size(Sarray));
EGCTPRtots = zeros(size(Sarray));
GCTPRs = zeros(size(Sarray));
PPVs = zeros(size(Sarray));
PPVtots = zeros(size(Sarray));
N = sum(EGC_seq_sizes{2});
M = sum(GC_seq_sizes{2});
for i = 1:length(Sarray)
    S = Sarray(i);

    EGCTPRtots(i) = TPRsameTtot(EGC_seq_sizes{2}, GC_seq_sizes{2}, GC_seq_sizes{1}, EGC_seq_sizes_neighbors, GC2_EGCmatch_seq_size_neighbors, S);
    GCTPRs(i) = TPRsameT(GC_seq_sizes{2}, EGC_seq_sizes{2}, GC_seq_sizes_neighbors{2}, S);
    EGCTPR2s(i) = TPRdiffT(EGC_seq_sizes{2}, GC2_EGCmatch_seq_size, GC2_EGCmatch_seq_size_neighbors, S);



%     PPVs(i) = (EGCTPRs(i)*N/(N+M))   /   ( (EGCTPRs(i)*N/(N+M)) + (GCTPRs(i)*M/(N+M)));
%     EGCTPRtots(i) = EGCTPR2s(i)+EGCTPR2s(i)-EGCTPRs(i)*EGCTPR2s(i);
    PPVtots(i) = (EGCTPRtots(i)*N/(N+M))   /   ( (EGCTPRtots(i)*N/(N+M)) + (GCTPRs(i)*M/(N+M)));
end

%% Simulation
SarraySim = [30, 50, 70];
N = sum(EGC_seq_sizes{2});
M = sum(GC_seq_sizes{2});
sensitivityarr = zeros(size(SarraySim));
precisionarr = zeros(size(SarraySim));
for k=1:length(SarraySim)
S = SarraySim(k);
repeat = 10^3;
sensitivity = zeros(1,repeat);
precision = zeros(1,repeat);
TP = zeros(1,repeat);
FN = zeros(1,repeat);
FP = zeros(1,repeat);
for i=1:repeat
    sample = struct('GC',{[],[]},'EGC',{[],[]});
    sample(1).GC = randsample(GCProperties{1},S);
    s1 = binornd(S,M/(N+M));
    sample(2).GC = randsample(GCProperties{2},s1);
    sample(2).EGC = randsample(EGCProperties{2},S-s1);
    [sensitivity(i), precision(i), TP(i), FN(i), FP(i)] =...
        labelingStats(sample,EGC_sequences,EGC_neighbors,GC_sequences,GC_neighbors);
end
sensitivity = sensitivity(~isnan(precision));
precision = precision(~isnan(precision));
sensitivityarr(k) = mean(sensitivity);
precisionarr(k) = mean(precision);
end
%%
% Plot
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
figure
t = tiledlayout(1,1,'Padding','tight');
t.Units = 'centimeters';
t.OuterPosition = [3, 3, 7, 7];
ax = nexttile;
plot(Sarray, EGCTPRtots, 'LineWidth', 1);
hold on
plot(Sarray, PPVtots, 'LineWidth', 1);
ylim([-0.1,1])
scatter(SarraySim, sensitivityarr, 30, colors{1});
scatter(SarraySim, precisionarr, 30, colors{2});
leg = legend({'Sensitivity (Bayesian)', 'Precision (Bayesian)', ...
    'Sensitivity (Monte Carlo)', 'Precision (Monte Carlo)'} ...
    ,'location','SouthEast');
leg.ItemTokenSize = [10,3];
xlabel('Number of sequences sampled')
ylabel('Sensitivity and precision')
%%
f = gcf;
savefig(f, fullfile('..','figures',[mfilename,'.fig']));
exportgraphics(f,fullfile('..','figures',[mfilename,'.pdf']),'ContentType','vector',...
            'BackgroundColor','none');
%%
function [sensitivity, precision, TP, FN, FP] = labelingStats(sample,EGC_sequences,EGC_neighbors,GC_sequences,GC_neighbors)
TP = 0;
FN = 0;
FP = 0;
for i=1:length(sample(2).EGC)
    if sum(sample(2).EGC(i)==sample(2).EGC)>1 || ismember(sample(2).EGC(i), sample(1).GC) ...
            || numneighbor(sample(2).EGC(i), sample(2).EGC, EGC_sequences{2}, EGC_neighbors)>1 ...
            || numneighbor(sample(2).EGC(i), sample(1).GC, EGC_sequences{2}, EGC_neighbors)>0
        TP = TP + 1;
    else
        FN = FN + 1;
    end
end
for i=1:length(sample(2).GC)
    if sum(sample(2).GC(i)==sample(2).GC)>1 ...
            || numneighbor(sample(2).GC(i), sample(2).GC, GC_sequences{2}, GC_neighbors{2})>1
        FP = FP + 1;
    end
end
sensitivity = TP/(TP+FN);
precision = TP/(TP+FP);
end

function numberOfNeighbors = numneighbor(sample_seq, samples, all_sequences, neighbors)
    idx = find(sample_seq==all_sequences);
    neighbors_seq = all_sequences(full(neighbors(idx,:))==1);
    numberOfNeighbors = sum(arrayfun(@(x) sum(samples==x), neighbors_seq));
end

function nums1TPR = TPRdiffT(nums1, nums2, nums2_neighbors, Samplenum)
N = sum(nums1);
L = sum(nums2);
Qs = Qprime(L,nums2_neighbors,Samplenum);
nums1TPR = 1 - sum(nums1/N.*Qs);
end

function nums1TPR = TPRsameT(nums1, nums2, nums1_neighbors, Samplenum)
N = sum(nums1);
M = sum(nums2);
Qs = Q(N, M, nums1_neighbors, Samplenum);
nums1TPR = 1 - sum(nums1/N.*Qs);
end

function nums1TPR = TPRsameTtot(EGCnums, GC2nums, GC1nums, EGC_nums_neighbors, GC1_nums_neighbors, Samplenum)
N = sum(EGCnums);
M = sum(GC2nums);
L = sum(GC1nums);
Qs = Q(N, M, EGC_nums_neighbors, Samplenum);
Qprimes = Qprime(L, GC1_nums_neighbors, Samplenum);
nums1TPR = 1 - sum(EGCnums/N.*Qs.*Qprimes);
end

function vals = Qprime(L,nums,S)
n = max(nums);
val = 1;
vals = zeros(1,length(nums));
for i=1:(n)
    val = val * ((L-S+1-i)/(L+1-i));
    if ismember(i, nums)
        idx = nums==i;
        vals(idx) = val;
    end
end
end

function vals = Q(N,M,nums,S)
n = max(nums);
val = 1;
vals = zeros(1,length(nums));
for i=1:(n)
    val = val * ((N+M-S+2-i)/(N+M+1-i));
    if ismember(i, nums)
        idx = nums==i;
        vals(idx) = val;
    end
end
end


