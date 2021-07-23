
classdef Portfolio
    properties
        finalAssetIdxList; %list of indexes of assets included in the portfolio
        initialAssetList; %assets listed in S&P500
        sectorSet; %set of sectors
        numStocksPerSec; 
        idxSec; %cell containing index of each asset associated with sector s
        numStocksPermcapQ;
        idxQ; %same as idxSec, but for mcapQ
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
            obj.numStocksPermcapQ = zeros(1,5); %indexed to sectorSet numStocksPerSec(s) = |obj.sectorSet(s)|
            obj.idxQ = cell(1,5);
            %count assets in each quintile
            for q = 1:5
                idx_q = [];
                for i = 1:length(obj.initialAssetList)
                    if(obj.initialAssetList{1,i}.mCapQ == q)
                        obj.numStocksPermcapQ(q) = obj.numStocksPermcapQ(q) + 1;
                        idx_q = [idx_q i];
                    end
                end
                obj.idxQ{1,q} = idx_q;
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