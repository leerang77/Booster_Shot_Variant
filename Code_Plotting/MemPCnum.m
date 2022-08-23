function [pcnum, memnum] = MemPCnum(param, result)
% Number of EGC and GC-derived plasma and memory B cells
buildup3D = @(x) reshape(x, size(x,1), [], param.n_ep);
pcnum_EGC = buildup3D(squeeze(result.output.pcnumbytime(2:2:end,:,2)));
memnum_EGC = buildup3D(squeeze(result.output.memnumbytime(2:2:end,:,2)));
pcnum_GC = buildup3D(squeeze(result.output.pcnumbytime(1:2:end-1, :, 2)));
memnum_GC = buildup3D(squeeze(result.output.memnumbytime(1:2:end-1, :, 2)));
for ep=1:2
   pcnum{ep} = [mean(pcnum_EGC(:,:,ep)); mean(pcnum_GC(:,:,ep))];
   memnum{ep} = [mean(memnum_EGC(:,:,ep)); mean(memnum_GC(:,:,ep))];   
end
end