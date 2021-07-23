%function [optPortfolios,RcumulOpt, overallCPUTimeOpt, benchRcumul, cplexobj] = ...
function [optPortfolios,RcumulOpt, overallCPUTimeOpt, benchRcumul] = ...
    getOptimalResults(model, Portfolios, benchPortfolios, lambda, period)
    T0 = cputime;
    optPortfolios = cell(1,period);
    %period = 2;
    cplexObj = cell(1,period);
    for t = 1:period
        assets = Portfolios{1,t}.initialAssetList;
        auxIdx = 1; %optportIDX
        if(t>1) %if it's not the first period 
            %start optimization
            %clc
            t0 = cputime;
            fprintf('\noptimising %s portfolio',Portfolios{1,t}.date);
            %%%%%%%%%SOLVE WITH GRASP-QUAD%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %solve for classical MVmodel
            if(strcmp('classical',model) )
                fprintf('(using standard MV) ...')
                %{
                [Portfolios{1,t}.w, Portfolios{1,t}.finalAssetIdxList, cplexobj{1, t-1}] = ...
                    CPLEX_classical(Portfolios{1,t}, lambda);
                %}
                [Portfolios{1,t}.w, Portfolios{1,t}.finalAssetIdxList] = ...
                   CPLEX_classical(Portfolios{1,t}, lambda);
                    %fprintf('\n|K| = %d', length(Portfolios{1,t}.finalAssetIdxList));
                %}
            end
            if(strcmp('TE',model) )
                fprintf('(using TE) ...')
                [Portfolios{1,t}.w, Portfolios{1,t}.finalAssetIdxList] = ...
                    CPLEX_TE(Portfolios{1,t}, lambda); 
            end
            %{
            %new graspquad
            [Portfolios{1,t}.w, ~,...
                Portfolios{1,t}.finalAssetIdxList] = core_graspQ(...
                Portfolios, optPortfolios, t, lambda, 20); 
            %}
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            optimizationCpuTime = cputime - t0;
            fprintf('\ndone - optimisation CPU time: %fs\n',optimizationCpuTime);
            %get optPortfolio
            Portfolios{1,t}.finalAssetIdxList
            Assets_opt_t = cell(1,length(Portfolios{1,t}.finalAssetIdxList));
            counti = 0;
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
    overallCPUTimeOpt = cputime - T0;
    %%Calculate performance metrics
    PerformanceMetrics = calcPerformance(optPortfolios, benchPortfolios, period);
    RcumulOpt = PerformanceMetrics{2,6};
    %benchRcumul = PerformanceMetrics{3,6};
    benchRcumul = PerformanceMetrics{2,3} %turnover
    
    %%Analyse Metrics
    clf
    x = 2:period;
    %{
    %Compare rbench_t and ropt_t
    hold on
    plot(x, PerformanceMetrics{3,2}, 'b');
    plot(x, PerformanceMetrics{2,2}, 'g');
    legend('rbench_IBOV', 'roptGRASP');
    title('Comparisson between roptGRASP and rbench_IBOV for all periods')
    %}
    %Compare K = 5, K= 10, K = 20 turnover
    hold on
    %plot(x, PerformanceMetrics{2,3}, 'b'); %turnover
    %plot(x, PerformanceMetrics{2,2}, 'g');
    %{  
    clf
    x = 1:period-1;
    hold on
    plot(x, turnover5, 'b'); %turnover
    plot(x, turnover10, 'g'); %turnover
    plot(x, turnover20, 'r'); %turnover
    legend('K = 5', 'K = 10', 'K = 20');
    xlabel('periodo')
    ylabel('turnover')
    %}
    %legend('rbench_IBOV', 'roptGRASP');
    %title('Comparisson between roptGRASP and rbench_IBOV for all periods')
end