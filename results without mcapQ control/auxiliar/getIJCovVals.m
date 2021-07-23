%%This function reads the covMat raw string
%%This funcition returns i j vectors of asset indexes and its associated cov values
%% co
function [i, j, cov] = getIJCovVals(rawCovMatFile, P)%search cov pairs and add it to y
i=zeros(1,length(rawCovMatFile{1,1})-1); %asset i index
j=zeros(1,length(rawCovMatFile{1,1})-1); %asset j index
cov = zeros(1,length(rawCovMatFile{1,1})-1); %vector of cov values
currentSEDOL = '';
current_i = 0;
lastpercentage = 0;
%{
Now some variables will be added to accelerate the reading process

How covMat file is arranged?
if we have 5 assets (with SEDOLS: A B C D E) in date YYYY-MM-DD, then the pattern is
YYYY-MM-DD A A 
YYYY-MM-DD A B
YYYY-MM-DD A C
YYYY-MM-DD A D
YYYY-MM-DD A E
YYYY-MM-DD B B
YYYY-MM-DD B C
YYYY-MM-DD B D
YYYY-MM-DD B E
YYYY-MM-DD C C
YYYY-MM-DD C D
YYYY-MM-DD C E
YYYY-MM-DD D D
YYYY-MM-DD D E
YYYY-MM-DD E E

we only need to search for index associated with SEDOLS in the first run
(while reading covs in relation to asset with SEDOL A)
%}
patternVector = zeros(1, length(P.initialAssetList)); %pattern vector containing index order of covMat file
auxIdx_lower = 1; %patternVector lower index 
auxIdx_variable = auxIdx_lower; %patternVecotr variable index
counter = 0;
for line=1:(length(rawCovMatFile{1,1})-1)
    rawLine = textscan(rawCovMatFile{1,1}{line+1}, '%s %s %s %s',... 
    'delimiter', ',');%obs.: line + 1 because the first line of the file is the header
    %rawLine{1,2} - SEDOL of i
    %rawLine{1,3} - SEDOL of j
    %rawLine{1,4} - covariance of i and j
    if(~strcmp(currentSEDOL, rawLine{1,2})) %compare currentSedol with i's SEDOL
        %fprintf('change i');
        %counter = 0 %jump file Header
        %counter = 1 %reading pattern
        %counter = 2 %pattern memorized
        counter = counter + 1;
        if(counter == 0)
            current_i = getIdx(rawLine{1,2}, P.initialAssetList)
        end
        if(counter < 2) %pattern not memorized
            current_i = getIdx(rawLine{1,2}, P.initialAssetList);
        else %index pattern memorized
            auxIdx_lower = auxIdx_lower + 1;
            auxIdx_variable = auxIdx_lower; %reset variable index
            current_i = patternVector(auxIdx_lower);
        end
        currentSEDOL = rawLine{1,2};
        %if the compared SEDOLS are not equal, this means we need
        %to change the current i index (asset changed)
    end
        if(counter < 2) %pattern not memorized
            i(line) = current_i;
            j(line) = getIdx(rawLine{1,3}, P.initialAssetList);
            patternVector(line) = j(line); %memorizing
        else %index pattern memorized
            i(line) = current_i;
            j(line) = patternVector(auxIdx_variable);
            auxIdx_variable = auxIdx_variable + 1;
        end
        cov(line) = str2double(rawLine{1,4});
        nowloading = (line/(length(rawCovMatFile{1,1})-1));
        if(nowloading > (lastpercentage*1.1))
            clc
            fprintf('loading covMat...%.2f%%', nowloading*100);
            lastpercentage = nowloading;
        end
end
end
