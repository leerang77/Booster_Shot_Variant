function [result,dirname] = load_result(p)
%% Define file name for saving the data
param = initializeParameters(p.vaxnum, p.T, p.k, p.numshot, p.E1h, p.dE12, p.p, ...
 p.masking, p.C0, p.w1, p.w2, p.steric, p.memToGCFrac,...
 p.outputprob, p.outputpcfrac, p.rho, p.earlybooster, p.tmax, p.first, p.last);

data_locs = {'Data_Prime', 'Data_Secondary', 'Data_Tertiary', 'Data_Vax4'};
for i=1:param.vaxnum
    if param.memToGCFrac > 0 && i>1
        data_locs{i} = fullfile('memToGC',sprintf('memToGCFrac_%.3f', ...
                                        param.memToGCFrac),data_locs{i});
    end
    if param.earlybooster && i==3
        data_locs{i} = fullfile('earlyBooster',data_locs{i});
    end
    if param.steric ~= 0
        data_locs{i} = fullfile('steric',sprintf('steric_%.2f', ...
                                        param.steric), data_locs{i});
    end
    if param.w1 > 0
        data_locs{i} = fullfile('agCaptureSaturation', data_locs{i});
    end
end
for i=1:param.vaxnum
    data_locs{i} = fullfile('Data', data_locs{i});
end
inputs = num2cell([param.vaxnum, param.E1h, param.dE12, param.p, ...
 param.masking, param.C0, param.w1, param.w2, param.steric,...
 param.outputprob, param.outputpcfrac, param.rho, param.tmax]);
inputs = cellfun(@(x) num2str(x), inputs, 'UniformOutput',false);
vaxtiming = [0, 28, 180, 180];

dirnm = cell(1,param.vaxnum);
fnm = cell(1,param.vaxnum);

[filepath, ~, ~] = fileparts(mfilename('fullpath'));
dirnm{param.vaxnum} = fullfile(filepath,'..',data_locs{param.vaxnum},...
                            strjoin(inputs(1:end),'_'));
timing = vaxtiming(param.vaxnum)*(param.earlybooster==0)+param.earlybooster;
if param.vaxnum>1
    dirnm{param.vaxnum-1} = fullfile('..',data_locs{param.vaxnum-1},...
     strjoin([num2str(param.vaxnum-1),inputs(2:end-1),num2str(timing)],'_'));
end
dirname = dirnm{param.vaxnum};

fprintf('\n%s',dirname)
if not(isfolder(dirname))
    error('Data Path Does Not Exist')
end
% if not(isfolder(result_dirname))
%     mkdir(result_dirname)
% end
result = loadResult(param.first, param.last, dirname);
result = combineResult(result, param.vaxnum);
end

function result = combineResult(resultarr, numDosing)
result.param = resultarr{1}.param;
result.conc = cellfun(@(obj) obj.conc, resultarr, 'UniformOutput', false);


affbytime = cellfun(@(obj) obj.gc.affbytime, resultarr, 'UniformOutput', false);
result.gc.affbytime = vertcat(affbytime{:});
numbytime = cellfun(@(obj) obj.gc.numbytime, resultarr, 'UniformOutput', false);
result.gc.numbytime = vertcat(numbytime{:});

naive = cellfun(@(obj) obj.naive, resultarr, 'UniformOutput', false);
result.naive = padvertcat(naive);

numbylineage = cellfun(@(obj) obj.gc.numbylineage, resultarr, 'UniformOutput', false);
result.gc.numbylineage = padvertcat(numbylineage);


finalpc = cellfun(@(obj) obj.output.finalpc, resultarr, 'UniformOutput', false);
result.output.finalpc = padvertcat(finalpc);
finalmem = cellfun(@(obj) obj.output.finalmem, resultarr, 'UniformOutput', false);
result.output.finalmem = padvertcat(finalmem);
finalpb = cellfun(@(obj) obj.output.finalpb, resultarr, 'UniformOutput', false);
result.output.finalpb = padvertcat(finalpb);


pcnumbytime = cellfun(@(obj) obj.output.pcnumbytime, resultarr, 'UniformOutput', false);
result.output.pcnumbytime = vertcat(pcnumbytime{:});
pcaffbytime = cellfun(@(obj) obj.output.pcaffbytime, resultarr, 'UniformOutput', false);
result.output.pcaffbytime = vertcat(pcaffbytime{:});
memnumbytime = cellfun(@(obj) obj.output.memnumbytime, resultarr, 'UniformOutput', false);
result.output.memnumbytime = vertcat(memnumbytime{:});
memaffbytime = cellfun(@(obj) obj.output.memaffbytime, resultarr, 'UniformOutput', false);
result.output.memaffbytime = vertcat(memaffbytime{:});
dead.plasmaCells = cellfun(@(obj) obj.dead.plasmaCells, resultarr, 'UniformOutput', false);
result.dead.plasmaCells = alignLeft(horzcat(dead.plasmaCells{:}), 1); 
if numDosing>1
    finalmemEGC = cellfun(@(obj) obj.memoryCellsEGC, resultarr, 'UniformOutput', false);
    result.finalmemEGC = alignLeft(horzcat(finalmemEGC{:}), 1);
    plasmaCellsEGC = cellfun(@(obj) obj.plasmaCellsEGC, resultarr, 'UniformOutput', false);
    result.plasmaCellsEGC = alignLeft(horzcat(plasmaCellsEGC{:}), 1);

    dead.plasmaCellsEGC = cellfun(@(obj) obj.dead.plasmaCellsEGC, resultarr, 'UniformOutput', false);
    result.dead.plasmaCellsEGC = alignLeft(horzcat(dead.plasmaCellsEGC{:}), 1);    
   
end
end

function result = loadResult(first, last, dirname)
% Combine smaller files from the folder the compile the data
% inputs:
%   inputs - parameters converted into cell array
%   first, last - range of the single GC run numbers
% outputs:
%   result(1xlast cell array): Each cell contains single GC run result

existing = zeros(1,last-first+1);
% dirname = fullfile('Data',dirname);
listing = dir(fullfile(dirname,'*.mat'));
listing = {listing.name};
result = cell(1,length(listing));
for i=1:length(listing)
   exp = '(?<first>\d*)_to_(?<last>\d*).mat';
   token = regexp(listing{i},exp,'names');
   n1 = round(str2double(token.first));
   n2 = round(str2double(token.last));
   if n1>=first && n2<=last
       existing(n1:n2) = 1;
       data = load(fullfile(dirname,listing{i}));
       result{i} = data.result;
       result{i}.first = n1;
   end
end
if any(existing==0)
   fprintf(num2str(find(existing==0)))
   error(['Missing files exist for: ', dirname]) 
end
nonempty = cellfun(@(x) ~isempty(x), result);
result = result(nonempty);
end