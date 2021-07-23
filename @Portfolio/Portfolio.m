
classdef Portfolio
    properties
        finalAssetIdxList; %list of indexes of assets included in the portfolio
        initialAssetList; %assets listed in S&P500
        sectorSet; %set of sectors
        numStocksPerSec; 
        idxSec; %cell containing index of each asset associated with sector s
        date; %time period associated with this portfolio
        solver; %solver that will optimize this portfolio
        %what we need to build the obj function? those 3 things:
        i_asset %index vector used to order cov index according to initialAssetList index order
        j_asset %index vector used to order cov index according to initialAssetList index order
        covij; %vector of covij values
        lambda; %risk aversion coefficient
        w %positions
        d %active weights
        %here you can add some performance metrics
    end
    methods
        function obj = Portfolio(Assets, d)
            obj.initialAssetList = Assets;
            obj.date = d;
            obj.i_asset = [];
            obj.j_asset = [];
            obj.covij = [];
            if(~isempty(obj.initialAssetList) > 0)
                %get sector set
                obj.sectorSet = obj.initialAssetList{1,1}.sector; %convert this array of chars to string
                for i = 2:length(obj.initialAssetList)
                   %compare strings!
                   currentSector = obj.initialAssetList{1,i}.sector;
                   if(~ismember(currentSector,obj.sectorSet))
                       obj.sectorSet = [obj.sectorSet currentSector]; %add new sector
                   end
                end
                obj.numStocksPerSec = zeros(1,length(obj.sectorSet)); %indexed to sectorSet numStocksPerSec(s) = |obj.sectorSet(s)|
                obj.idxSec = cell(1,length(obj.sectorSet));
                %count assets in each sector
                for s = 1:length(obj.sectorSet)
                    idx_s = [];
                    for i = 1:length(obj.initialAssetList)
                        currentSector = obj.initialAssetList{1,i}.sector;
                        if(ismember(currentSector,obj.sectorSet(s)))
                            obj.numStocksPerSec(s) = obj.numStocksPerSec(s) + 1;
                            idx_s = [idx_s i];
                        end
                    end
                    obj.idxSec{1,s} = idx_s; 
                end
            end
        end
    end
    methods(Static)
        solver = setSolver(solverName);
        postAssetList = optimizePortfolio(P, rawCovMatFile); %returns a list of assets
        %[cumulativeReturnP, cumulativeReturnB] = getCumRetP() %get
        %cumulative return of P and Benchmark
        %infoRatio = getIR() % get Information Ratio
        %here you can put other performance metric functions
    end
end