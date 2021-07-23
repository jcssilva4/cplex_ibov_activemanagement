function [sector_satisfied, mcapq_satisfied] = AdjustSecMcapQ(S, sec_ok, mcap_ok)
    %S must contain all sectors and mcapq 
    firstcheck = 1; %first time checking mcap and sec
    proceedthis = 0; %depends on sec_ok, mcap_ok
    while(~proceedthis)
        if(first) %first time running
            [sec_ok, mcap_ok] = checkSecMcapQ(S, sec_ok, mcap_ok); %check if each one is satisfied
            if(sec_ok + mcap_ok > 1) %check if both are satisfied
                proceedthis = 1; %S contains all sectors and mcapq 
                fprintf('S contains all sectors and mcapq');
            end
        end
        if(sec_ok + mcap_ok < 2) %check if both are satisfied
            s = randsample(length(RCL),1); %not satisfied, sample another stock
            Stemp = [S RCL(s)]; %add a new stock to try to change (msec_ok + mcap_ok)
            sectemp_ok = 0;
            mcapq_ok = 0;
            [sectemp_ok, mcapqtemp_ok] = checkSecMcapQ(Stemp, sec_ok, mcap_ok);
            if (sectemp_ok - sec_ok) %detect change
                S = Stemp; %add a new stock
                sec_ok = 1;
            end
             if (mcapqtemp_ok - mcapq_ok) %detect change
                S = Stemp; %add a new stock
                mcapq_ok = 1;
            end
        else  %both are satisfied
             proceedthis = 1; %S contains all sectors and mcapq 
             fprintf('S contains all sectors and mcapq');
        end
    end
    if(~sec_ok)
        %get sector set
        sectorSet = P.initialAssetList{1,S(1)}.sector; %convert this array of chars to string
        for i = 2:length(S)
           %compare strings!
           currentSector = P.initialAssetList{1,S(i)}.sector;
           if(~ismember(currentSector,sectorSet))
               sectorSet = [sectorSet currentSector]; %add new sector
           end
        end
        if(~(length(sectorSet) < length(P.sectorSet)))
            sector_satisfied = 1; %S contains all sectors
        end
    else
        sector_satisfied = 1; %S contains all sectors already
    end
    if(~mcap_ok)
    end
end