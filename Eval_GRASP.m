%% Evaluate GRASPs
%% READ DATA, SETUP SIMULATIONS
clear
clc
format long
addpath([pwd '/auxiliar/'])
addpath(['C:/Program Files/IBM/ILOG/CPLEX_Studio128/cplex/matlab/x64_win64/'])
addpath([pwd '/CPLEX/'])
addpath([pwd '/GRASP_QUAD/'])
addpath([pwd '/PerformanceMetrics/']);

% Read data
inputData = Data(); %create a Data object
%inputDataPath = 'C:\Users\GREEFO-CAA\Desktop\IbovespaData\';
inputDataPath = [pwd '\IbovespaData\'];
%read IBOVESPA timeseries.txt
inputData.tsSet = Data.readSP500file([inputDataPath 'timeseries\ibovespaTimeSeries.txt']); clc
%Backtesting setup 
auxIdx = 2; %auxiliar index (get assets name, dates, asset betas)
counter = 1;
Portfolios = cell(1,21); %initialize portfolio variable
nRebalancesPerYear = zeros(1,2); idxRebalance = 1;
initialPeriod = 1;
finalPeriod = 21;
period = initialPeriod;
while(period <= finalPeriod) 
    %get Ibovespa assets
    [initialAssetList, auxIdx, refDate] = Asset.getAll(inputData.tsSet, auxIdx);
    %read the covMat file associated with refDate
    fprintf ('reading covmatrix data from %s...', refDate); 
    %getcovMatrix of yyyy-mm-dd
    fid = fopen([inputDataPath 'covmat\covmat_', refDate, '.txt'], 'r');
    raw = textscan(fid, '%s %s %s', 'delimiter', '\t'); 
    fclose(fid);
    fprintf ('ok\n');
    %create a portfolio instance
    P = Portfolio(initialAssetList, refDate); %initialize a new portfolio P
    [P.i_asset, P.j_asset, P.covij] = getIJCovVals(raw, P); %this function is contained in \auxiliar\
    Portfolios{1,period} = P;
    period = period + 1;
end
%get initial portfolio (0 assets)
P = Portfolio(initialAssetList, 'initial'); %initial portfolio
[P.i_asset, P.j_asset, P.covij] = getIJCovVals(raw, P); %this function is contained in \auxiliar\
Portfolios{1,period} = P;
%flip portfolios vertically (begin from 2016-05-20)
tempP = Portfolios;
for p = 1:period
    Portfolios{1,p} = tempP{1,period-p+1};
end
%%get Bench portfolios
benchPortfolios = cell(1,period-1);
for t = 1:period
    Assets_bench_t = cell(1,length(Portfolios{1,t}.initialAssetList));
    assets = Portfolios{1,t}.initialAssetList;
    auxbIdx = 1; %benchportIDX
    if(t>1) %if it's not the first period 
        for i = 1:length(assets)
            %get benchPortoflio
            if(assets{1,i}.benchWeight > 0)
                Assets_bench_t{1,auxbIdx} = Asset(assets{1,i}.sedol,...
                    assets{1,i}.name, assets{1,i}.sector, assets{1,i}.benchWeight,...
                        assets{1,i}.alphaScore, assets{1,i}.r, assets{1,i}.benchWeight); %check the last arg
                auxbIdx = auxbIdx + 1;
            end
        end
    else %if it is the first period (base portfolio)
        for i = 1:length(assets)
            %get benchPortoflio
            if(assets{1,i}.benchWeight > 0)
                Assets_bench_t{1,auxbIdx} = Asset(assets{1,i}.sedol,...
                     assets{1,i}.name, assets{1,i}.sector, assets{1,i}.benchWeight,...
                        assets{1,i}.alphaScore, assets{1,i}.r, assets{1,i}.benchWeight); %check the last arg
                auxbIdx = auxbIdx + 1;
            end
        end
    end
    benchPortfolios{1,t} = Portfolio(Assets_bench_t, Portfolios{1,t}.date);
end

%% get Roptimal - classical
clc
lambda = 1;
%{
 [optPortfolios, RcumulCPLEX, overallCPUTimeCPLEX, benchRcumul cplex] = ...
        getOptimalResults('classical', Portfolios, benchPortfolios, lambda,...
        period);
%}

[optPortfolios,RcumulCPLEX, overallCPUTimeCPLEX, benchRcumul] = ...
        getOptimalResults('classical', Portfolios, benchPortfolios, lambda,...
        period);
    %}
    turnover5 = benchRcumul;
%% Plot turnover

%{
CplexCPUTime = zeros(1,length(cplex));
CplexGAP = zeros(1,length(cplex));
CplexObj = zeros(1,length(cplex));
for t = 1:length(cplex)
    CplexCPUTime(t) = cplex{1,t}.Solution.time;
    CplexGAP(t) = cplex{1,t}.Solution.miprelgap;
    CplexObj(t) = cplex{1,t}.Solution.bestobjval;
end
Results_CPLEX{1,1} = 'Gaps'; Results_CPLEX{1,2} = 'bestobjval'; 
Results_CPLEX{1,3} = 'MeanCPUTime(s)'; 
Results_CPLEX{2,1} = CplexGAP; Results_CPLEX{2,2} = CplexObj; 
Results_CPLEX{2,3} = mean(CplexCPUTime); 
        %}
%% get Roptimal - Tracking error
clc
lambda = 1;

 [optPortfolios, RcumulCPLEX, overallCPUTimeCPLEX, benchRcumul] = ...
        getOptimalResults('TE', Portfolios, benchPortfolios, lambda,...
        period);
%{
CplexCPUTime = zeros(1,length(cplex));
CplexGAP = zeros(1,length(cplex));
CplexObj = zeros(1,length(cplex));
for t = 1:length(cplex)
    CplexCPUTime(t) = cplex{1,t}.Solution.time;
    CplexGAP(t) = cplex{1,t}.Solution.miprelgap;
    CplexObj(t) = cplex{1,t}.Solution.bestobjval;
end
Results_CPLEX_TE{1,1} = 'Gaps'; Results_CPLEX_TE{1,2} = 'bestobjval'; 
Results_CPLEX_TE{1,3} = 'MeanCPUTime(s)'; 
Results_CPLEX_TE{2,1} = CplexGAP; Results_CPLEX_TE{2,2} = CplexObj; 
Results_CPLEX_TE{2,3} = mean(CplexCPUTime); 
        %}
