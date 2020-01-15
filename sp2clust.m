function clust = sp2clust(sp, ksDir)

clear clust;
nClu = numel(sp.clusterId);
for iS = 1:nClu
    clust(iS) = clustClass(sp, sp.clusterId(iS));
end

disp('Done!')


%% save clust
if exist('ksDir', 'var')
    save(fullfile(ksDir, 'clust.mat'), 'clust')
    disp('Done saving ''clust''')
else
    disp('no ''ksDir'' provided as input so I''m not saving nada')
end

