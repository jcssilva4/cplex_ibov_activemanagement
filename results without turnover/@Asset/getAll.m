function [list, eIdx] = getAll(tsSet, returns, bIdx)
    %bIdx -> initial index position
    %eIdx -> final index position
    listIdx = 1; %assists asset list construction
    auxIdx = bIdx; %assists asset count
    referenceDate = tsSet{1,1}{bIdx};
    %count the number of assets composing S&P500 for that period
    while(strcmp(referenceDate, tsSet{1,1}{auxIdx}))
        auxIdx = auxIdx + 1;
        if(auxIdx>length(tsSet{1,1})) %avoids to exceed matrix dimension
            break;
        end
    end
    eIdx = auxIdx; %get final index position
    list = cell(1, eIdx - bIdx); %initialize asset list 
    for i = bIdx:(eIdx-1)
        list{1, listIdx} = Asset(tsSet{1,2}{i}, tsSet{1,3}{i},...
            str2double(tsSet{1,4}{i}),str2double(tsSet{1,5}{i}),tsSet{1,6}{i}...
            , str2double(tsSet{1,7}{i}), str2double(tsSet{1,8}{i}), returns(listIdx),0);
        listIdx = listIdx + 1;
    end
end