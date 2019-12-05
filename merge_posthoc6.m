function rez = merge_posthoc6(rez)

fpThresh = 0.01;
ops = rez.ops;
Nfilt = rez.ops.Nfilt;

tfi = rez.iNeigh;
tf = rez.cProj;

clusterIDs = rez.st3(:,2);

nSpikes = size(rez.st3,1);

fmax = zeros(nSpikes,1, 'single');
ppairs = cell(Nfilt,1); % potential pairs
fpRate   = ones(Nfilt,1);
nspikes = zeros(Nfilt,1);
depths  = rez.ypos(:);  
refDur = 30;
temps = reshape(rez.dWU, [], Nfilt); 

for testID = 1:Nfilt
    spikesTest = clusterIDs==testID;
    
    pp = tfi(:, testID);
    pp(pp==testID) = [];
    ppairs{testID} = pp;
    [~, isame] = min( abs(tfi(:, testID)-testID));
    fmax(spikesTest, 1) = tf(spikesTest, isame);
    
    nspikes(testID) = sum(spikesTest);
    if nspikes(testID)==0
        continue
    end
    fpRate(testID) = ISIViolations(rez.st3(spikesTest,1), 0, refDur);
end


%% Step 1: only merge clusters with no fp violations
% clear iMegaC
iMegaC = {};
new_spikes = [];
picked = zeros(Nfilt,1);
picked(fpRate > fpThresh) = true; % eliminate false positive violations

inewclust = 0;

valid = find(~picked);

[~, order] = sort(nspikes(valid), 'descend');

clustList = valid(order);


for iseed = 1:numel(clustList)
    
    fprintf('New Cluster Total: %d\n', inewclust)
    this = clustList(iseed);
    
    if picked(this)
        continue
    end
    
    run_list = this;
    
    do_not_consider = [];
    
    all_possible = setdiff(tfi(:, find(sum(tfi==this))), this);
    all_possible = union(ppairs{this}, all_possible(:));
    all_possible = intersect(all_possible, valid);
    if exist('iMegaC', 'var')
        all_possible = setdiff(all_possible, unique(cell2mat(iMegaC)));
    end
    
    pair_list = all_possible(abs(depths(all_possible)-depths(this)) < 1);
    consider = pair_list;
    
    %     consider = pair_list(fpRate(pair_list) < fpThresh);
%     do_not_consider = [do_not_consider; setdiff(all_possible, consider)];
    
    strun = find(clusterIDs==this);
    
    % npairs to consider merging
    nc = numel(consider);
    
    %% first to consider
    for iconsider = 1:nc
        ipair = consider(iconsider);
        
        imm = ismember(tfi(:, ipair), run_list);
        
        if nspikes(ipair) < 100 % don't merge new units that don't spike
            do_not_consider = [do_not_consider; ipair];
            continue
        end
        
        rho = corr(temps(:, [run_list ipair]));
        rho = mean(rho(~triu(rho)));
        if rho < .5
            do_not_consider = [do_not_consider; ipair];
            continue
        end
        
        if ~any(imm)
            do_not_consider = [do_not_consider; ipair];
            continue
        end
        
        
        
        
        new_spikes = find(clusterIDs==ipair);
        f1new = max(tf(new_spikes, imm), [], 2);
        f2new = fmax(new_spikes);
        
        f1old = fmax(strun);
        
        f2old = NaN * ones(numel(f1old), 1, 'single');
        
        if ~ismember(tfi(:, run_list), ipair)
            feats = intersect(tfi(:,ipair), tfi(:,run_list));
        else
            feats = ipair;
        end
        
        
        i0 = 0;
        for j = 1:length(run_list)
            ifeat = find(ismember(tfi(:, run_list(j)), setdiff(feats, run_list(j))),1);
            if ~isempty(ifeat)
                ix = clusterIDs == run_list(j);
                f2old(i0 + (1:nspikes(run_list(j)))) = tf(ix, ifeat);
                i0 = i0 + nspikes(run_list(j));
            end
        end
        
        f1old(isnan(f2old))=[];
        f2old(isnan(f2old))=[];
        
        figure(2); clf
        plot(f1old, f2old, '.')
        hold on
        plot(f1new, f2new, '.')
        drawnow
        
        if isempty(f1old)
            do_not_consider = [do_not_consider; ipair];
            continue
        end
        
        mo = merging_score(f1old - f2old, f1new-f2new, ops.fracse);
        
        if mo > 3
            do_not_consider = [do_not_consider; ipair];
            continue
        end
        
%         fprintf('Checking xcorr\n')
%         % for xcorr analysis
%         binSize = 30;
%         nlags = 100;
%         s1 = rez.st3(strun,1);
%         s1 = ceil(s1/binSize);
%         
%         
%         s2 = rez.st3(new_spikes,1);
%         s2 = ceil(s2/binSize);
%         
%         num = max(max(s1), max(s2));
%         
%         sp1 = full(sparse(s1, ones(numel(s1),1), ones(numel(s1),1), num, 1));
%         sp2 = full(sparse(s2, ones(numel(s2),1), ones(numel(s2),1), num, 1));
%         
%         xc3 = xcorr(sp1, sp2, nlags, 'unbiased');
%         if sum(sp1) > sum(sp2)
%             xc1 = xcorr(sp1, nlags, 'unbiased');
%         else
%             xc1 = xcorr(sp2, nlags, 'unbiased');
%         end
%         
%         xc1(nlags+1) = 0;
%         xc1 = smooth(xc1, 5);
%         xc3 = smooth(xc3, 5);
%                 
%         nfun = @(x) x/norm(x);
%         xc1 = nfun(xc1);
%         xc3 = nfun(xc3);
%         
%         figure(1); clf
%         lags = -nlags:nlags;
%         plot(lags,xc1, lags, xc3)
%         xproj = corr([xc1 xc3]);
%         title(xproj(2))
%         drawnow
%         
%         
%         if xproj(2) > .7
            
            fprintf('Merging %d\n', ipair)
            
            strun = cat(1, strun, new_spikes);
            run_list(end+1) = ipair;
            picked(ipair)   = 1;
            
%             pair_list = unique(cat(1, pair_list, ppairs{ipair}));
%             pair_list(ismember(pair_list, run_list)) = [];
            
%         end
        
    end
    
    if numel(run_list) > 1
        inewclust = inewclust + 1;
        iMegaC{inewclust} = run_list;
    end
end

%% now merge high FP rate units
picked = zeros(Nfilt,1);
picked(unique(cell2mat(iMegaC))) = 1;
picked(fpRate < .2) = 1; % eliminate low fp clusters

valid = find(~picked);

[~, order] = sort(nspikes(valid), 'descend');

clustList = valid(order);


for iseed = 1:numel(clustList)
    
    fprintf('New Cluster Total: %d\n', inewclust)
    this = clustList(iseed);
    
    if picked(this)
        continue
    end
    
    run_list = this;
    
    do_not_consider = [];
    
    all_possible = setdiff(tfi(:, find(sum(tfi==this))), this);
    all_possible = union(ppairs{this}, all_possible(:));
    all_possible = intersect(all_possible, valid);
    all_possible = setdiff(all_possible, unique(cell2mat(iMegaC)));
    
    pair_list = all_possible(abs(depths(all_possible)-depths(this)) < 1);
    consider = pair_list;
    
%     strun = find(clusterIDs==this);
    
    % npairs to consider merging
    nc = numel(consider);
    
    %% first to consider
    for iconsider = 1:nc
        
        strun = find(ismember(clusterIDs, run_list));
        
        ipair = consider(iconsider);
        
        imm = ismember(tfi(:, ipair), run_list);
                
        rho = corr(temps(:, [run_list ipair]));
        rho = mean(rho(~triu(rho)));
        
        if rho > .5
            fprintf('Merging %d\n', ipair)
            
            strun = cat(1, strun, new_spikes);
            run_list(end+1) = ipair;
            picked(ipair)   = 1;
            
            pair_list = unique(cat(1, pair_list, ppairs{ipair}));
            pair_list(ismember(pair_list, run_list)) = [];
            continue
            
        end
            
        if ~any(imm)
            do_not_consider = [do_not_consider; ipair];
            continue
        end
        
        
        
        
        new_spikes = find(clusterIDs==ipair);
        f1new = max(tf(new_spikes, imm), [], 2);
        f2new = fmax(new_spikes);
        
        f1old = fmax(strun);
        
        f2old = NaN * ones(numel(f1old), 1, 'single');
        
        if ~ismember(tfi(:, run_list), ipair)
            feats = intersect(tfi(:,ipair), tfi(:,run_list));
        else
            feats = ipair;
        end
        
        
        i0 = 0;
        for j = 1:length(run_list)
            ifeat = find(ismember(tfi(:, run_list(j)), setdiff(feats, run_list(j))),1);
            if ~isempty(ifeat)
                ix = clusterIDs == run_list(j);
                f2old(i0 + (1:nspikes(run_list(j)))) = tf(ix, ifeat);
                i0 = i0 + nspikes(run_list(j));
            end
        end
        f2old = f2old(1:numel(f1old));
        
        f1old(isnan(f2old))=[];
        f2old(isnan(f2old))=[];
        
        figure(2); clf
        plot(f1old, f2old, '.')
        hold on
        plot(f1new, f2new, '.')
        drawnow
        
        if isempty(f1old)
            do_not_consider = [do_not_consider; ipair];
            continue
        end
        
        mo = merging_score(f1old - f2old, f1new-f2new, ops.fracse);
        
        if mo > 3
            do_not_consider = [do_not_consider; ipair];
            continue
        end
        
%         fprintf('Checking xcorr\n')
%         % for xcorr analysis
%         binSize = 30;
%         nlags = 100;
%         s1 = rez.st3(strun,1);
%         s1 = ceil(s1/binSize);
%         
%         
%         s2 = rez.st3(new_spikes,1);
%         s2 = ceil(s2/binSize);
%         
%         num = max(max(s1), max(s2));
%         
%         sp1 = full(sparse(s1, ones(numel(s1),1), ones(numel(s1),1), num, 1));
%         sp2 = full(sparse(s2, ones(numel(s2),1), ones(numel(s2),1), num, 1));
%         
%         xc3 = xcorr(sp1, sp2, nlags, 'unbiased');
%         if sum(sp1) > sum(sp2)
%             xc1 = xcorr(sp1, nlags, 'unbiased');
%         else
%             xc1 = xcorr(sp2, nlags, 'unbiased');
%         end
%         
%         xc1(nlags+1) = 0;
%         xc1 = smooth(xc1, 5);
%         xc3 = smooth(xc3, 5);
%                 
%         nfun = @(x) x/norm(x);
%         xc1 = nfun(xc1);
%         xc3 = nfun(xc3);
%         
%         figure(1); clf
%         lags = -nlags:nlags;
%         plot(lags,xc1, lags, xc3)
%         xproj = corr([xc1 xc3]);
%         title(xproj(2))
%         drawnow
%         
%         
%         if xproj(2) > .7
            
            fprintf('Merging %d\n', ipair)
            
            strun = cat(1, strun, new_spikes);
            run_list(end+1) = ipair;
            picked(ipair)   = 1;
            
            pair_list = unique(cat(1, pair_list, ppairs{ipair}));
            pair_list(ismember(pair_list, run_list)) = [];
            
%         end
        
    end
    
    if numel(run_list) > 1
        inewclust = inewclust + 1;
        iMegaC{inewclust} = run_list;
    end
end

%%
% grouped = unique(cell2mat(iMegaC));
% notgrouped = setdiff(1:Nfilt, grouped);
% 
% nMC = numel(iMegaC);
% notgroupedid = (1:numel(notgrouped)) + nMC;

fprintf('Removing Collisions from merge\n')
noiseCluster = iMegaC{end}(end);

refDur = 30;
iMega = zeros(Nfilt, 1);
for i = 1:length(iMegaC)
    fprintf('%d/%d\n', i, length(iMegaC))
    iMega(iMegaC{i}) = iMegaC{i}(1);
    
    pairs = nchoosek(iMegaC{i},2);
    if all(fpRate(iMegaC{i}) < fpThresh)
        removeIx = [];
        for j = 1:size(pairs,1)
            % check for collisions upon merge
            % only take one spike when both spiked at the same time
            ix1 = find(rez.st3(:,2)==pairs(j,1));
            ix2 = find(rez.st3(:,2)==pairs(j,2));
            t1 = rez.st3(ix1,1);
            t2 = rez.st3(ix2,1);
            
            n1 = numel(t1);
            n2 = numel(t2);
            removeix = [];
            if n1 < n2
                for jt = 1:n1
                    removeix = [removeix; ix2(abs(t1(jt) - t2) < refDur)];
                    %                     t2(abs(t1(jt) - t2) < refDur) = inf;
                end
            else
                for jt = 1:n2
                    removeix = [removeix; ix1(abs(t2(jt) - t1) < refDur)];
                end
            end
        end
        
        % asign collisions to noise Cluster
        rez.st(removeix,2) = noiseCluster;
    end
end

rez.iMega = iMega;
rez.iMegaC = iMegaC;


rez.st3(:,5) = iMega(rez.st3(:,2));