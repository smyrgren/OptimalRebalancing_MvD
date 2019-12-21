clear classes

%% Create the MvD object
mvd = OptReb_MvD();

% Initial - Portfolio and asset characterstics
astNames = ["AA" "BB" "CC"];
curWts = [.4 .4 .2];
benWts = [.6 .3 .1];
expRts = [.16 .06 .11];
stdDvs = [.19 .08 .12];
corMtx = [1 -.6 -.2; -.6 1 .5; -.2 .5 1];
covMtx = corMtx.*(stdDvs'*stdDvs);
tnsCst = 1 * [.0050 .0030 .0080];

% Cost function
funcName = 'MVT';
lambdaAbs = 10;
lambdaRel = 0;

% Map of properties
propMap = containers.Map({'assetNames', 'curWts', 'benWts', 'expRts', 'covMtx', 'tnsCst', ...
                         'funcName', 'lambdaAbs', 'lambdaRel'}, ...
                         {["AA" "BB" "CC"], curWts, benWts, expRts, covMtx, tnsCst, ...
                           funcName, lambdaAbs, lambdaRel});
mvd = setProperties(mvd, propMap);         

%% Regular Optimization - Step 1
% Optimal portfolio
mvd = optimizeCostfunction(mvd, 1);

% As close as we can get with transaction costs
mvd = optimizeCostfunction(mvd, 0);
                       
%% MvD - Approximation of best multiperiod rebalancing solution
simFreq = 12;
numPeriods = 36;
numSims = 100;


tic
fprintf('Training/Calibration  - MvD - Quadratic Approximation ... \n')
for iMvD = 1:7
    
    mvd(end+1) = mvd(1);
    
    QAprx = .1 * (iMvD - 1) * eye(length(curWts));

    % Map of properties
    propMap = containers.Map({'funcName', 'QAprx', 'simFreq', 'numPeriods', 'numSims'}, ...
                             {'QAprx', QAprx, simFreq, numPeriods, numSims});
    mvd(end) = setProperties(mvd(end), propMap);
    
    mvd(end) = trainingSimulation(mvd(end), iMvD + 1);

end
toc

% ´The best single-period heuristic to approximate a multi-period optimization
% problem
f1 = figure;
hold
xlabel('QAPrx #');
ylabel('Cost (bps)')
plot([1:length(mvd)-1], [mvd(2:end).avgCEC] + [mvd(2:end).avgTC], 'ro')
bar([1:length(mvd)-1], [mvd(2:end).avgTC], 'b');
legend('Total Cost (bps)', 'Transaction costs (bps)');
hold off


%% Mvd Approximation - Simulation
totCost = [mvd(2:end).avgCEC] + [mvd(2:end).avgTC];
indxBest = find(totCost == min(totCost));
bestMvD = mvd(1 + indxBest);

% We need a greater number of simulations
numSims = 500;
propMap = containers.Map({'numSims'}, {numSims});
bestMvD = setProperties(bestMvD, propMap);

% Run the rebalance simulation
bestMvD = rebalanceSimulation(bestMvD, [1 3 6 12 40], .01 * [1 2 3 4 5]);

% Extract the comparative stats
statsTab = extractStatistics(bestMvD);

