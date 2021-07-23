%% Algorithm used for Ibovespa Backtesting

%% READ DATA, SETUP SIMULATIONS
clear
clc
format shortg
addpath([pwd '/auxiliar/'])
addpath(['C:/Program Files/IBM/ILOG/CPLEX_Studio128/cplex/matlab/x64_win64/'])
addpath([pwd '/CPLEX/'])
addpath([pwd '/GRASP_QUAD/'])
addpath([pwd '/PerformanceMetrics/']);

% Read data
inputData = Data(); %create a Data object
inputDataPath = 'C:\Users\GREEFO-CAA\Desktop\IbovespaData\';
%inputDataPath = 'C:\Users\Julio\Desktop\IbovespaData\';
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
%% Optimise Portfolios
clc
T0 = cputime;
optPortfolios = cell(1,period);
lambda = 1;
for t = 1:period
    assets = Portfolios{1,t}.initialAssetList;
    auxIdx = 1; %optportIDX
    if(t>1) %if it's not the first period 
        %start optimization
        %clc
        t0 = cputime;
        fprintf('\noptimising %s portfolio...',Portfolios{1,t}.date);
        %%%%%%%%%SOLVE WITH GRASP-QUAD%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %old graspquad
        [~, Portfolios{1,t}.w,...
            ~,Portfolios{1,t}.finalAssetIdxList] = core_graspQ(...
            Portfolios{1,t}, lambda, 100); 
        %{
        %new graspquad
        [Portfolios{1,t}.w, ~,...
            Portfolios{1,t}.finalAssetIdxList] = core_graspQ(...
            Portfolios, optPortfolios, t, lambda, 20); 
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        optimizationCpuTime = cputime - t0;
        fprintf('done - optimisation CPU time: %fs',optimizationCpuTime);
        %get optPortfolio
        Assets_opt_t = cell(1,length(Portfolios{1,t}.finalAssetIdxList));
        for i = 1:length(assets)    
            if(ismember(i,Portfolios{1,t}.finalAssetIdxList))
                Assets_opt_t{1,auxIdx} = Asset(assets{1,i}.sedol,...
                    assets{1,i}.name, assets{1,i}.sector, assets{1,i}.benchWeight,...
                    assets{1,i}.alphaScore, assets{1,i}.r, Portfolios{1,t}.w(i));
                auxIdx = auxIdx + 1;
            end
        end
    else %if it is the first period (initial portfolio)
        %Assets_opt_t = [];
        Portfolios{1,t}.w = zeros(1,length(assets));
        Assets_opt_t = cell(1,length(Portfolios{1,t}.initialAssetList));
        for i = 1:length(assets)
            if(Portfolios{1,t}.w(i) >= 0) %check this condition...only wi is required
                Assets_opt_t{1,auxIdx} = Asset(assets{1,i}.sedol,...
                    assets{1,i}.name, assets{1,i}.sector, assets{1,i}.benchWeight,...
                    assets{1,i}.alphaScore, assets{1,i}.r, 0);
                auxIdx = auxIdx + 1;
            end
        end
    end
    optPortfolios{1,t} = Portfolio(Assets_opt_t, Portfolios{1,t}.date); %initialize a new portfolio P
end
overallCpuTime = cputime - T0

%% get Bench portfolios
clc
benchPortfolios = cell(1,period-1);
for t = 1:period
    Assets_bench_t = cell(1,length(Portfolios{1,t}.initialAssetList));
    assets = Portfolios{1,t}.initialAssetList;
    auxbIdx = 1; %benchportIDX
    if(t>1) %if it's not the first period 
        for i = 1:length(assets)
            %get benchPortoflio
            Assets_bench_t{1,auxbIdx} = Asset(assets{1,i}.sedol,...
                assets{1,i}.name, assets{1,i}.sector, assets{1,i}.benchWeight,...
                    assets{1,i}.alphaScore, assets{1,i}.r, assets{1,i}.benchWeight); %check the last arg
            auxbIdx = auxbIdx + 1;
        end
    else %if it is the first period (base portfolio)
        for i = 1:length(assets)
            %get benchPortoflio
            Assets_bench_t{1,auxbIdx} = Asset(assets{1,i}.sedol,...
                 assets{1,i}.name, assets{1,i}.sector, assets{1,i}.benchWeight,...
                    assets{1,i}.alphaScore, assets{1,i}.r, assets{1,i}.benchWeight); %check the last arg
            auxbIdx = auxbIdx + 1;
        end
    end
    benchPortfolios{1,t} = Portfolio(Assets_bench_t, Portfolios{1,t}.date);
end
%% Calculate performance metrics
clc
PerformanceMetrics = calcPerformance(optPortfolios, benchPortfolios, period);
    %PerformanceMetrics = calcPerformance(optPortfolios, benchPortfolios, 124+1);
%% Write Results
%HERE

%% Analyse Metrics
clf
x = 2:period;
%Compare rbench_t and ropt_t
hold on
plot(x, PerformanceMetrics{3,2}, 'b');
plot(x, PerformanceMetrics{2,2}, 'g');
legend('rbench_IBOV', 'roptGRASP');
title('Comparisson between roptGRASP and rbench_IBOV for all periods')
%{
%% 
clf
%Compare turnoverbench_t and turnoveropt_t
%plot(x, PerformanceMetrics{3,3}, 'b');
hold on
plot(x, PerformanceMetricsNOTO{2,3}, 'k');
plot(x, PerformanceMetricsTO{2,3}, 'r');
legend('TOoptNOTO', 'TOoptWITHTO');
title('Comparisson of TO for all periods')
%% 
clf
%compare cumulRbench and cumulRopt
hold on
plot(x, PerformanceMetricsNOTO{2,3}, 'k');
plot(x, PerformanceMetricsTO{2,3}, 'r');
legend('TOoptNOTO', 'TOoptWITHTO');
title('Comparisson of TO for all periods')
%}