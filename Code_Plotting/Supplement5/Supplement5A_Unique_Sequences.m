fprintf(mfilename)
addpath(fullfile('..','..','Code_Parameter_Generation'));
addpath(fullfile('..','..','Code_Result_Analysis'));
addpath(fullfile('..','..','Code_Simulation'));
addpath('..');
%%
% GC
p = base_case_parameters();
p.last = 200;
[result{1},~] = load_result(p);
p.vaxnum = 2; p.tmax = 180;
[result{2},~] = load_result(p);
GC_seq_sizes = cell(1,2); EGC_seq_sizes = cell(1,2); 
GCProperties = cell(1,2); EGCProperties = cell(1,2);
f = 7;
neighborMaxDist = 1;
%%
GC_neighbors = cell(1,2); GC_sequences = cell(1,2);
for i=1:2
    GC_neighbors{i} = sparse(10^5,10^5);
    output = result{i}.output.finalmem;
    n_GC = size(output,1);
    num_uniques = [];
    GCProperties{i} = zeros(1,10^6);
    GC_sequences{i} = [];
    idx = 1;
    seqidx = 1;
    for k=1:n_GC
       mem = squeeze(output(k,output(k,:,3)>f,:));
       if ~iscolumn(mem)
           mem = mem';
       end
       if ~isempty(mem)
       Property = sum(mem([1,3,4,9],:));
       GCProperties{i}(idx:idx+length(Property)-1) = Property;
       idx = idx + length(Property);
       [sequences, uniqueIdx] = unique(Property);
       nums = histcounts(Property, [sequences,inf]);
       GC_sequences{i} = [GC_sequences{i}, sequences];
       num_uniques = [num_uniques, nums];
       neighbors = getNeighbors(mem, uniqueIdx,neighborMaxDist);
       GC_neighbors{i}(seqidx:seqidx+length(sequences)-1,...
           seqidx:seqidx+length(sequences)-1) = neighbors;
       seqidx = seqidx + length(sequences);
       end
    end
    GCProperties{i} = GCProperties{i}(GCProperties{i}>0);
    [GC_seq_sizes{i},Idx] = sort(num_uniques,'descend');
    GC_neighbors{i} = GC_neighbors{i}(Idx,Idx);
    GC_sequences{i} = GC_sequences{i}(Idx);
    GC_seq_sizes_neighbors{i} = zeros(size(GC_seq_sizes{i}));
    for j=1:length(GC_neighbors{i})
        GC_seq_sizes_neighbors{i}(j) = sum(GC_seq_sizes{i}(GC_neighbors{i}(j,:)==1));
    end
end

%% EGC
mem = result{2}.finalmemEGC;
Property = sum(mem([1,3,4,9],:));
Property = Property(mem(3,:)>f);
[sequences, uniqueIdx] = unique(Property);
% EGC_neighbors = getNeighbors(mem, uniqueIdx, 1);
cnts = arrayfun(@(x) sum(x==Property), sequences);
[EGC_seq_sizes{2},Idx] = sort(cnts, 'descend');
EGC_sequences{2} = sequences(Idx);

id = arrayfun(@(x) find(GC_sequences{1}==x), EGC_sequences{2});
EGC_neighbors = GC_neighbors{1}(id,id);

GC2_EGCmatch_seq_size = arrayfun(@(x) sum(x==GCProperties{1}), sequences);
GC2_EGCmatch_seq_size = [GC2_EGCmatch_seq_size(Idx)];

EGCProperties{2} = Property;
EGC_seq_sizes_neighbors = zeros(size(EGC_seq_sizes{2}));
GC2_EGCmatch_seq_size_neighbors = zeros(size(GC2_EGCmatch_seq_size));
for j=1:length(EGC_neighbors)
    EGC_seq_sizes_neighbors(j) = sum(EGC_seq_sizes{2}(EGC_neighbors(j,:)==1));
    GC2_EGCmatch_seq_size_neighbors(j) = sum(GC2_EGCmatch_seq_size(EGC_neighbors(j,:)==1));
end
outputFile = fullfile('..','Results',[mfilename,'_',num2str(neighborMaxDist),'.mat']);
save(outputFile, 'GC_seq_sizes', 'EGC_seq_sizes', 'EGC_seq_sizes_neighbors',...
    'GC_seq_sizes_neighbors', 'GC2_EGCmatch_seq_size', 'GC2_EGCmatch_seq_size_neighbors', 'GCProperties',...
    'EGCProperties', 'GC_neighbors', 'EGC_neighbors', 'GC_sequences', 'EGC_sequences')

function neighbors = getNeighbors(mem, uniqueIdx, maxdist)
n = length(uniqueIdx);
lineage = transpose(squeeze(mem(1,uniqueIdx)));
criteria = cell(1,2);
criteria{1} = zeros(n);
for i=1:2
    mut{i} = int2bit(squeeze(mem(7+i-1,uniqueIdx)),40);
    for j=1:n
        for k=j+1:n
            criteria{1}(j,k) = ...
                criteria{1}(j,k) + sum(xor(mut{i}(:,j),mut{i}(:,k)));
            criteria{1}(k,j) = criteria{1}(j,k);
        end
    end
end
criteria{1} = criteria{1} <= maxdist;
criteria{2} = repmat(lineage,1,n) == transpose(repmat(lineage,1,n));
neighbors = criteria{1} & criteria{2};
end