function A = padvertcat(A)
%A is a cell array
lengths = cellfun(@(x) size(x,2), A);
maxlength = max(lengths);
if any(lengths~=maxlength)
    A = cellfun(@(x) [x, zeros(size(x,1),maxlength-size(x,2),size(x,3))], A, 'UniformOutput', false);
    A = vertcat(A{:});
else
    A = vertcat(A{:});
end
end