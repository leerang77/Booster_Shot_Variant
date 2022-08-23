function createParameters(txt_name, writeoption, p, varargin)
% Inputs: 
%   txt_name - Name of the parameter text file with or without '.txt'
%                 'a' for appending to the existing content

if not(isfolder(fullfile('..','parameters')))
    mkdir(fullfile('..','parameters'))
end

txt_name = erase(txt_name, '.txt');
fid = fopen([txt_name,'.txt'], writeoption);

% Write the text file for simulation runs
first_arr = linspace(p.first, p.last+1, p.numfrag+1);
first_arr(end)=[];
numrepeat = (p.last-p.first+1)/p.numfrag;

for first = first_arr
    last = first+numrepeat-1;
    if isempty(varargin)
    inputs = num2cell([p.vaxnum, p.T, p.k, p.numshot, p.E1h, p.dE12, p.p, ...
          p.masking, p.C0, p.w1, p.w2, p.steric, p.memToGCFrac,p.outputprob, ...
          p.outputpcfrac, p.rho, p.earlybooster, p.tmax, first, last]);
    else
    inputs = num2cell([p.T, p.k, p.numshot, p.E1h, p.dE12, p.p, ...
          p.masking, p.C0, p.w1, p.w2, p.steric, p.memToGCFrac,p.outputprob, ...
          p.outputpcfrac, p.rho, p.earlybooster, first, last]);
    end
    inputs = cellfun(@(x) num2str(x), inputs, 'UniformOutput', false);
    fprintf(fid, [strjoin(inputs(1:end),'\t'),'\n']);
end
fclose(fid);
end