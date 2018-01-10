
nSmall = 10;

%%
idxGood     = false(numel(tSamples,1)); % init
for ii = 1:nSmall
    ptrUp(ii) = find(tSamples >= tUp(ii),1);
    ptrDn(ii) = find(tSamples >= tDn(ii),1);
    idxGood(ptrUp(ii):ptrDn(ii)) = true;
end
tt = tSamples(idxGood);

%%

zz = [];
for ii = 1:nSmall
    
    zz = [zz tSamples(sUp(ii):sDn(ii))];
    
end


%%
figure, hold on

plot(tt)
plot(zz)