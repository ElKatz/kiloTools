function clust = sp2clust(sp, ksDir)


nClu = numel(sp.clusterId);
for iS = 1:nClu
    clusterId = sp.clusterId(iS);
    clust(iS) = clustClass;
    clust(iS).sp2clust_forClusterId(sp, clusterId);
end

disp('Done!')


%% save clust
if exist('ksDir', 'var')
    save(fullfile(ksDir, 'clust.mat'), 'obj')
    disp('Done saving ''clust''')
else
    disp('no ''ksDir'' provided as input so I''m not saving nada')
end

