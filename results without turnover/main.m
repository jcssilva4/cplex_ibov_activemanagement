%% Algorithm that solves 2018 INFORMS OR & Analytics Student Team Competition problem

%% READ DATA, SETUP SIMULATIONS
clear
clc
format shortg
addpath([pwd '/auxiliar/'])
%addpath([pwd '/MILXPM/'])
%addpath([pwd '/MILXPMsimple/'])
addpath(['C:/Program Files/IBM/ILOG/CPLEX_Studio128/cplex/matlab/x64_win64/'])
addpath([pwd '/CPLEX/'])
%addpath('C:/Program Files/ampl_pack/matbar/Interface')
%addpath([pwd '/MIDACO Solver/'])
addpath([pwd '/GRASP_QUAD/'])
addpath([pwd '/PerformanceMetrics/']);

dDay = 7*4; %time variation between simulations (in days)
currDay = 3;
currMonth = 1;
currYear = 2007;
lastYear = 2007;
%{
finalDay = 28;
finalMonth = 03;
finalYear = 2007;
%}

finalDay = 21;
finalMonth = 12;
finalYear = 2016;
%}
first = 1; %first time executing?
stop = 0;
auxIdx = 2; %auxiliar index (get assets name, dates, asset betas)
counter = 1;
Portfolios = cell(1,132); %initialize portfolio variable
nRebalancesPerYear = zeros(1,10); idxRebalance = 1;
% Read data
inputData = Data(); %create a Data object
%inputDataPath = 'C:\Users\GREEFO-CAA\Desktop\sp500\inputData\';
inputDataPath = 'C:\Users\Julio\Desktop\matlab_prototype\inputData\';
%read S&P500 timeseries file
inputData.tsSet = Data.readSP500file([inputDataPath 'Timeseries_data_SP500.txt']); clc
%read Results Template file (get historical returns)
inputData.returns = Data.readResultsTempFile([inputDataPath 'returns.csv']); clc
% Start periodic optimization
%{
    remember that S&P500 is an index composed by ~500 stocks
    This composition changes from period to period (stocks listed
    are not constant)
%}
period = 2;
%Simulations setup -> starting period: 31-01-2007
while (~stop) 
    [initialAssetList, auxIdx] = Asset.getAll(inputData.tsSet,...
        inputData.returns, auxIdx);%get S&P500 assets
    if(~first) %change this condition so that the algorithm can start in different periods(detect changes in date)
        [currYear, currMonth, currDay] = getNextDate(lastYear, lastMonth, lastDay, dDay);
    else
        %get initial portfolio
        P1 = Portfolio(initialAssetList, 'baseportfolio'); %initialize a new portfolio P
        P1.w = Data.readInitPortFile([inputDataPath 'initialPortfolio.txt']);
        Portfolios{1,1} = P1;
        first = 0; 
    end
    %read the covMat file associated with this date
    date = formatDate(currYear, currMonth, currDay);
    fprintf ('reading covmatrix data from %s...', date); 
    %getcovMatrix of yyyy-mm-dd
    fid = fopen([inputDataPath 'cov_mat_' date '.csv'], 'rt');
    raw = textscan(fid, '%s'); 
    fclose(fid);
    fprintf ('ok\n');
    %create a portfolio instance
    P = Portfolio(initialAssetList, sprintf('%d/%d/%d', currMonth,...
        currDay, currYear)); %initialize a new portfolio P
    %P.solver = Portfolio.setSolver('GA');
    %P.finalAssetList = Portfolio.optimizePortfolio(P, raw); %optimize P
    [P.i_asset, P.j_asset, P.covij] = getIJCovVals(raw, P);
    Portfolios{1,period} = P;
    if(lastYear~=currYear)%detect change in year
        idxRebalance = idxRebalance + 1;
    end
    nRebalancesPerYear(idxRebalance) = nRebalancesPerYear(idxRebalance) + 1; 
    period =  period + 1;
    %checking stop conditions...
    if(currDay == finalDay)
        if(currMonth == finalMonth)
            if(currYear == finalYear)
                stop = 1; 
            end
        end
    end
    lastYear = currYear;
    lastMonth = currMonth;
    lastDay = currDay;
    clc
end
%% Optimize Portfolios  
T0 = cputime;
for t = 2:period-1
    t0 = cputime;
    clc
    fprintf('\noptimizing %s portfolio\n',Portfolios{1,t}.date);
    %%%%%%%%%SOLVE WITH GRASP-QUAD%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Portfolios{1,t}.d, Portfolios{1,t}.w, ~,...
        Portfolios{1,t}.finalAssetIdxList] = core_graspQ(Portfolios{1,t}, 1, 400); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%GRASP PERFORMANCE%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %checkGQPerformance([10 50 100]) %plot compare UEF and CCEF generated
    %by GRASP, this function is included in matlab_prototy/auxiliar!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %get w's contained in x
    %get sigmas contained in x
    optimizationCpuTime = cputime - t0
end
overallCpuTime = cputime - T0
%% get Optimal and Bench portfolios
clc
optPortfolios = cell(1,period-1);
benchPortfolios = cell(1,period-1);
for t = 1:period-1
    Assets_opt_t = cell(1,length(Portfolios{1,t}.finalAssetIdxList));
    Assets_bench_t = cell(1,length(Portfolios{1,t}.initialAssetList));
    assets = Portfolios{1,t}.initialAssetList;
    auxIdx = 1; %optportIDX
    auxbIdx = 1; %benchportIDX
    if(t>1) %if it's not the first period (03-01-2007)
        %get optPortfolio
        for i = 1:length(assets)
            if(ismember(i,Portfolios{1,t}.finalAssetIdxList))
                Assets_opt_t{1,auxIdx} = Asset(assets{1,i}.sedol,...
                    assets{1,i}.sector, assets{1,i}.beta, assets{1,i}.alphaScore,...
                    assets{1,i}.name, assets{1,i}.benchWeight,...
                    assets{1,i}.mCapQ, assets{1,i}.r,Portfolios{1,t}.w(i));
                auxIdx = auxIdx + 1;
            end
            %get benchPortoflio
            Assets_bench_t{1,auxbIdx} = Asset(assets{1,i}.sedol,...
                assets{1,i}.sector, assets{1,i}.beta, assets{1,i}.alphaScore,...
                assets{1,i}.name, assets{1,i}.benchWeight,...
                assets{1,i}.mCapQ, assets{1,i}.r, assets{1,i}.benchWeight); %check the last arg
            auxbIdx = auxbIdx + 1;
        end
    else %if it is the first period (base portfolio)
        for i = 1:length(assets)
            if(Portfolios{1,t}.w(i) > 0) %check this condition...only wi is required
                Assets_opt_t{1,auxIdx} = Asset(assets{1,i}.sedol,...
                    assets{1,i}.sector, assets{1,i}.beta, assets{1,i}.alphaScore,...
                    assets{1,i}.name, assets{1,i}.benchWeight,...
                    assets{1,i}.mCapQ, assets{1,i}.r, Portfolios{1,t}.w(i));
                auxIdx = auxIdx + 1;
            end
            %get benchPortoflio
            Assets_bench_t{1,auxbIdx} = Asset(assets{1,i}.sedol,...
                assets{1,i}.sector, assets{1,i}.beta, assets{1,i}.alphaScore,...
                assets{1,i}.name, assets{1,i}.benchWeight,...
                assets{1,i}.mCapQ, assets{1,i}.r, assets{1,i}.benchWeight); %check the last arg
            auxbIdx = auxbIdx + 1;
        end
    end
    benchPortfolios{1,t} = Portfolio(Assets_bench_t, Portfolios{1,t}.date);
    optPortfolios{1,t} = Portfolio(Assets_opt_t, Portfolios{1,t}.date); %initialize a new portfolio P
end
%% Calculate performance metrics
clc
PerformanceMetrics = calcPerformance(optPortfolios, benchPortfolios, period);
%PerformanceMetrics = calcPerformance(optPortfolios, benchPortfolios, 3);
%% Write Results
%HERE

%% Analyse Metrics
clf
x = 1:period-2;
%Compare rbench_t and ropt_t
hold on
plot(x, PerformanceMetrics{3,2}, 'b');
plot(x, PerformanceMetrics{2,2}, 'g');
legend('rbench', 'roptGRASP');
title('Comparisson between roptGRASP and rbench for all periods')
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