function strobedEvents = getStrobedFromRaw(rawPath)

        % we'll get:
        %   tsMap - time vector
        %   evTs - event timestamp (i.e. strobes)
        %   evSv - event strobe value
        
        
% read the strobed word info (values & time stamps).
strobedEvents.eventInfo = PL2EventTs(rawPath, 'Strobed');
strobedEvents.RSTARTInfo = PL2EventTs(rawPath, 'RSTART');
strobedEvents.RSTOPInfo = PL2EventTs(rawPath, 'RSTOP');
