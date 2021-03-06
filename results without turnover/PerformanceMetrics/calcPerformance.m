function PerformanceMetrics = calcPerformance(optPortfolios, benchPortfolios, period)
nMetrics_t = 3; %metrics that vary with time (ropt, turnover, radjusted)
nMetrics = nMetrics_t + 6;
PerformanceMetrics = cell(3,nMetrics+1);
PerformanceMetrics{1,1} = 'Metrics';
PerformanceMetrics{2,1} = 'portfolio';
PerformanceMetrics{3,1} = 'benchmark';
PerformanceMetrics{1,2} = 'ropt_t';
PerformanceMetrics{1,3} = 'turnover_t';
PerformanceMetrics{1,4} = 'rTxad_t';
PerformanceMetrics{1,5} = 'adjustedInfoRatio'; PerformanceMetrics{2,5} = 0;
PerformanceMetrics{1,6} = 'Rcumul'; PerformanceMetrics{2,6} = 0;
PerformanceMetrics{1,7} = 'Ranual'; PerformanceMetrics{2,7} = 0;
PerformanceMetrics{1,8} = 'Rexcess_anual'; PerformanceMetrics{2,8} = 0;
PerformanceMetrics{1,9} = 'annualizedTE'; PerformanceMetrics{2,9} = 0;
PerformanceMetrics{1,10} = 'SharpeRatio'; PerformanceMetrics{2,10} = 0;
%calculate ropt + rbench + turnover + radjusted
for t = 2:nMetrics_t+1
    PerformanceMetrics{2,t} = zeros(1,period-2); %PerformanceMetrics{2,i} -  portfolio performance metric i for period t
    PerformanceMetrics{3,t} = zeros(1,period-2); %PerformanceMetrics{3,i} -  benchmark performance metric i for period t
end
IR_pt1_num = 1; IR_pt2_num = 1; %first and second parts of the adjIR numerator
for t = 2:period-1
    fprintf('\nperiod %d: computing performance metrics...\n',t-1)
    fprintf('computing overall return...');
    [PerformanceMetrics{2,2}(t-1),PerformanceMetrics{3,2}(t-1)] = get_ropt_t(optPortfolios{1,t},benchPortfolios{1,t});
    fprintf('ok!\ncomputing turnover...');
    PerformanceMetrics{2,3}(t-1) = get_turnover_t(optPortfolios{1,t-1}, optPortfolios{1,t});
    PerformanceMetrics{3,3}(t-1) = get_turnover_t(benchPortfolios{1,t-1}, benchPortfolios{1,t}); %bench
    fprintf('ok!\ncomputing adjusted return...');
    PerformanceMetrics{2,4}(t-1) = PerformanceMetrics{2,2}(t-1) - 0.005*PerformanceMetrics{2,3}(t-1);
    PerformanceMetrics{3,4}(t-1) = PerformanceMetrics{3,2}(t-1) - 0.005*PerformanceMetrics{3,3}(t-1); %bench
    fprintf('ok!\ncomputing fractions of the numerator of the adj information ratio...');
    %IR_pt1_num = IR_pt1_num*(1+PerformanceMetrics{2,4}(t-1));%*(1+radj_t)
    %IR_pt2_num = IR_pt2_num*(1+PerformanceMetrics{3,4}(t-1));%*(1+radjbench_t)
    IR_pt1_num = IR_pt1_num*(1+PerformanceMetrics{2,2}(t-1));%*(1+ropt_t)
    IR_pt2_num = IR_pt2_num*(1+PerformanceMetrics{3,2}(t-1));%*(1+rbench_t)
    fprintf('ok!\n');
end
std_RadjustedT_minus_RbenchT = std(PerformanceMetrics{2,4} - PerformanceMetrics{3,2});
%getIR
PerformanceMetrics{2,5} = (IR_pt1_num - IR_pt2_num)/std_RadjustedT_minus_RbenchT;
%get Rcumul_portfolio + Rcumul_bench 
PerformanceMetrics{2,6} = IR_pt1_num-1; %port
PerformanceMetrics{3,6} = IR_pt2_num-1; %bench
%get Ranual_port + Ranual_bench (## CHECK NUM OF REBALANCES PER YEAR!!!!!!!!!!!)
PerformanceMetrics{2,7} = PerformanceMetrics{2,6}; %port
PerformanceMetrics{3,7} = PerformanceMetrics{3,6}; %bench
%get Rexcess_anual 
PerformanceMetrics{2,8} = PerformanceMetrics{2,7} - PerformanceMetrics{3,7}; %port
%get annualizedTE
PerformanceMetrics{2,9} = sqrt(13)*std_RadjustedT_minus_RbenchT;
%get Sharpe Ratio
PerformanceMetrics{2,10} = PerformanceMetrics{2,6}/std(PerformanceMetrics{2,4});
end