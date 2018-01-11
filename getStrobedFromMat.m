function strobed = getStrobedFromMat(matFileName)
%% load stim codes list
disp('Loading stimcodes_FST');
codes = stimcodes_FST;

%% load Events from mat file
disp('Loading event ts & sv');
load(matFileName, 'evTs', 'evSv');
ts = evTs;
sv = evSv;

%% get strobed
disp('extracting trial event codes');
% strobed = get_strobed(sv,codes);
start = find(sv == codes.trialBegin);
stop =  find(sv == codes.trialEnd);
j = 1; strobed = {};
if (length(start) ~= length(stop)) || any((stop-start)<=0)
    disp('mismatched start and stop codes but carrying on fixing');
end
while ~isempty(start) && ~isempty(stop)
    if length(start) >=2 && length(stop) >= 1
        if start(1) < stop(1) && stop(1) < start(2)
            strobed{j} = [sv(start(1):stop(1)) ts(start(1):stop(1))];
            j = j+1;
            start(1)=[];stop(1)=[];
        elseif stop(1) < start(1) % start code missed
            stop(1) = [];
        elseif start(2) < stop(1) % stop code missed
            start(1) = [];
        end
    elseif length(start) == 1 && length(stop) == 1 && start(1) < stop(1)
        strobed{j} = [sv(start(1):stop(1)) ts(start(1):stop(1))];
        j = j+1;
        start(1)=[];stop(1)=[];
    else
        break
    end
end
