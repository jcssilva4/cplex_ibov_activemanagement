function [S, Sc, proceedGRASP] = adjustSecMcapQ(P, S, Sc, RCL)
    S = S;
    Sc = Sc; 
    proceedGRASP = 0;
    %S must contain all sectors and mcapq 
    [sec_ok, mcap_ok, sectorSet] = checkSecMcapQ(P, S, 0, 0); %get logical conditions
    if(sec_ok + mcap_ok > 1) %check if both are satisfied
        proceedGRASP = 1; %S contains all sectors and mcapq 
        %fprintf('S contains all sectors and mcapq');
    else
        %fprintf('Trynig to obtain all sectors and mcapq...');
        if(~sec_ok) %search the missing sectors in RCL
            for s=1:length(RCL)
                currentSector = P.initialAssetList{1,RCL(s)}.sector;
                if(~ismember(currentSector,sectorSet))
                   sectorSet = [sectorSet currentSector]; %add new sector
                   S = [S RCL(s)]; %add a new stock containing one of the remaining sectors
                   Sc = setdiff(Sc, RCL(s)); %Sc - s

                end
            end
        end
        [~, mcap_ok, ~, mcapqSet] = checkSecMcapQ(P, S, sec_ok, mcap_ok); %maybe we've found some missing mcapq while finding missing sectors
    end
    if(sec_ok + mcap_ok > 1) %check if both are satisfied
        proceedGRASP = 1;
        %fprintf('S contains all sectors and mcapq');
    else
        %fprintf('S does not contain all sectors and mcapq');
    end
end

function [sector_satisfied, mcapq_satisfied, sectorSet, mcapqSet] = checkSecMcapQ(P, S, sec_ok, mcap_ok)
    sector_satisfied = 0;
    mcapq_satisfied = 0;
    sectorSet = [];
    mcapqSet = [];
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
        sector_satisfied = 1; %S already contains all sectors
    end
end