function [Rcumul, overallCPUTime, benchRcumul, f] = ...
    getResults(model, Portfolios, benchPortfolios, lambda, iters, period)
    T0 = cputime;
    optPortfolios = cell(1,period);
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
                [~, Portfolios{1,t}.w,...
                    f,Portfolios{1,t}.finalAssetIdxList] = GQ_classical(...
                    Portfolios{1,t}, lambda, iters); 
                %if you change your model (no sector -> using sector)
                %you need to modify solve_classical(include ineq linear constraints)
                %and GRASP GQ classical_nosec -> GQ_classical
            end
            if(strcmp('classical_nosec',model) )
                fprintf('(using standard MV) ...')
                [~, Portfolios{1,t}.w,...
                    f,Portfolios{1,t}.finalAssetIdxList] = GQ_classical_no_sec(...
                    Portfolios{1,t}, lambda, iters); 
                %if you change your model (no sector -> using sector)
                %you need to modify solve_classical(include ineq linear constraints)
                %and GRASP GQ classical_nosec -> GQ_classical
            end
            if(strcmp('TE',model) )
                fprintf('(using TE) ...')
                [~, Portfolios{1,t}.w,...
                    f,Portfolios{1,t}.finalAssetIdxList] = GQ_TE(...
                    Portfolios{1,t}, lambda, iters); 
            end
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
    overallCPUTime = cputime - T0;

    %%Calculate performance metrics
    PerformanceMetrics = calcPerformance(optPortfolios, benchPortfolios, period);
    Rcumul = PerformanceMetrics{2,6};
    benchRcumul = PerformanceMetrics{3,6};
end