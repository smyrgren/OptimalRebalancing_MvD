classdef OptReb_MvD < Costfunction
    properties
        QAprx
        numSims
        numPeriods
        simFreq
        trnSim
        avgTC
        avgCEC
        strNames
        rebSim
        rebStats
    end
    methods (Access = public)
        
        % Constructor
        function obj = OptReb_MvD()
            obj@Costfunction();
        end
        
        % Set class properties
        function obj = setProperties(obj, propMap)
            
            objProps = properties(obj);
            mapKeys = propMap.keys;
            for ikey = 1:length(mapKeys)
                
                % Set the properties given  
                if any(strcmp(mapKeys(ikey), objProps))
                    obj.(mapKeys{ikey}) = propMap(mapKeys{ikey});
                else
                    fprintf('"%s" - is not an object property\n', mapKeys{ikey})
                end
   
            end
            
        end % setProperties
        
        function obj = trainingSimulation(obj, nMvD)
            
            % Reste the simulation structure
            obj.trnSim = [];

            % Portfolio input parameters
            expRts = (1 + obj.expRts).^(1 / obj.simFreq) - 1;
            covMtx = obj.covMtx / obj.simFreq;
            
            % Reset the rnadom number generator
            % rng(1)
            
            for is = 1:obj.numSims
                
                % Input weights
                iniWts = obj.optWts;
                
                fprintf('Session: %s \t Simulation#%s of %s\n', ...
                               int2str(nMvD - 1), int2str(is), int2str(obj.numSims))
            
                % Run through the periods 
                for ip = 1:obj.numPeriods
                    
                    % Generate the periodic returns
                    rndRts = mvnrnd(expRts, covMtx);

                    % Update the weights based on the returns over the period
                    endWts(ip,:) = (1 + rndRts).*iniWts(ip,:);
                    endWts(ip,:) = endWts(ip,:) / sum(endWts(ip,:));

                    % Update the portfolio object
                    propMap = containers.Map({'curWts'},{endWts(ip,:)});
                    obj = setProperties(obj, propMap);

                    % Optimize the weights basd on the csotfunction and the
                    % quadratic approximation
                    obj = optimizeCostfunction(obj, 0);

                    % The resulting optimal weights are the new initial weights
                    % in the next period
                    iniWts(ip+1,:) = obj.newWts;

                    % Turnover and transaction costs
                    turnOver(ip,1) = sum(abs(iniWts(ip+1,:) - endWts(ip,:)));
                    transCost(ip,1) = obj.tnsCst * abs(iniWts(ip+1,:) - endWts(ip,:))';
                    % crtEqDiff(ip,1) = obj.funcValCur- obj.funcValOpt;
                    crtEqDiff(ip,1) = obj.funcValNew - obj.funcValOpt; 

                end % numPeriods
                
                % Save the weights
                %obj.trnSim(is).iniWts = iniWts(2:end,:);
                %obj.trnSim(is).endWts = endWts;

                % Overall stats for the simulation
                obj.trnSim(is).turnOver = sum(turnOver);
                obj.trnSim(is).transCost = 1e4 * sum(transCost);
                obj.trnSim(is).crtEqCost = 1e4 * sum(crtEqDiff); 
                obj.trnSim(is).allPositive = all(crtEqDiff > 10*eps);
                
                
            end % numSims
            
           obj.avgTC = mean([obj.trnSim.transCost]);
           obj.avgCEC = mean([obj.trnSim.crtEqCost]);
           
        end % trainingSimulation
        
        % Simlulate rebalance ewith frequency and tolerance bands
        function obj = rebalanceSimulation(obj, rebFreqs, rebTols)
            
            % Reset the simualtion substruct
            obj.rebSim = [];
            obj.rebStats = [];

            % Portfolio input parameters
            expRts = (1 + obj.expRts).^(1 / obj.simFreq) - 1;
            covMtx = obj.covMtx / obj.simFreq;
            
            % How many portfolios are we running
            numFreqs = length(rebFreqs);
            numTols = length(rebTols);
            numPrt = 1 + numFreqs + numTols;
            
            % Reset the random number generator
            % rng(1)
            
            % Names rebalance stratgies
            colNames(1) = "QAprx";
            for ip = 1:numFreqs
                if rebFreqs(ip) > obj.numPeriods
                    colNames(1 + ip) = "No Rebalance";
                else
                    colNames(1 + ip) = "Freq@ - " + int2str(rebFreqs(ip));
                end
            end
            for ip = 1:numTols
                colNames(1 + numFreqs + ip) = "Tol@ - " + int2str(100 * rebTols(ip)) + "%";
            end
            obj.strNames = colNames;
            
            fprintf('Rebalance simulation ...\n')
            for is = 1:obj.numSims
                
                % Set initial variables
                iniWts = [];
                endWts = [];
                for ip = 1:numPrt
                    iniWts(1,:,ip) = obj.optWts;
                    obj.rebSim(is,ip).aum_nTC = 100;
                    obj.rebSim(is,ip).aum_wTC = 100;
                end
                
                % Progress?
                fprintf('Simulation#%s of %s\n', int2str(is), int2str(obj.numSims))

                % Run through the periods 
                % Really ought not save the individual simulations
                for period = 1:obj.numPeriods

                    % Generate the periodic returns
                    rndRts = mvnrnd(expRts, covMtx);

                    % Update the weights based on the returns over the period
                    for ip = 1:numPrt
                        endWts(period,:, ip) = (1 + rndRts).*iniWts(period,:, ip);
                        endWts(period,:, ip) = endWts(period,:,ip) / sum(endWts(period,:, ip));
                    end
                    
                    % Update the portfolio 
                    propMap = containers.Map({'curWts'},{endWts(period,:,1)});
                    obj = setProperties(obj, propMap);

                    % Optimize the weights based on the costfunction and the
                    % MvD quadratic approximation
                    obj = optimizeCostfunction(obj, 0);
                    
                    % Determine new weights
                    % The resulting optimal weights are the new initial weights
                    % in the next period
                    iniWts(period+1,:, 1) = obj.newWts;
                    
                    % Frequency rebalance
                    for ix = 1:numFreqs
                        if mod(period, rebFreqs(ix)) == 0
                            iniWts(period+1,:,1+ix) = obj.optWts;
                        else
                            iniWts(period+1,:,1+ix) = endWts(period,:,1+ix);
                        end
                    end
                    
                    % Tolerance bands
                    for ix = 1:numTols
                        if any(abs(endWts(period,:,1+numFreqs+ix) - obj.optWts) > rebTols(ix))
                            % Add half rebalance for comparison ... ?
                            iniWts(period+1,:,1+numFreqs+ix) = obj.optWts;
                        else
                            iniWts(period+1,:,1+numFreqs+ix) = endWts(period,:,1+numFreqs+ix);
                        end
                    end
                    
                    % Tunrover, transaction costs and aum
                    for ip = 1:numPrt
                        obj.rebSim(is, ip).turnOver(period,1) = sum(abs(iniWts(period+1,:,ip) - endWts(period,:,ip)));
                        obj.rebSim(is, ip).transCost(period,1) = obj.tnsCst*abs(iniWts(period+1,:,ip) - endWts(period,:,ip))';
                        obj.rebSim(is, ip).crtEqCost(period,1) = cfMVT(iniWts(period+1,:,ip)', obj.benWts', obj.expRts, ...
                                                                    obj.lambdaAbs, obj.lambdaRel, obj.covMtx) ...
                                                                    - obj.funcValOpt;
                        obj.rebSim(is, ip).aum_nTC(period+1,1) = (1 + iniWts(period,:,ip)*rndRts')*obj.rebSim(is, ip).aum_nTC(period,1);
                        obj.rebSim(is, ip).aum_wTC(period+1,1) = (1 + iniWts(period,:,ip)*rndRts' - ...
                                                                      obj.rebSim(is, ip).transCost(period,1))*obj.rebSim(is, ip).aum_wTC(period,1);
                    end
            
                end % end numPeriods
                
                for ip = 1:numPrt
                   obj.rebStats.TRC(is, ip) = 1e4 * sum([obj.rebSim(is,ip).transCost],1);
                   obj.rebStats.CEC(is, ip) = 1e4 * sum([obj.rebSim(is,ip).crtEqCost],1);
                   obj.rebStats.TOT(is, ip) = obj.rebStats.TRC(is, ip) + obj.rebStats.CEC(is, ip);
                   obj.rebStats.aum_nTC(is, ip) = obj.rebSim(is, ip).aum_nTC(end,1);
                   obj.rebStats.aum_wTC(is, ip) = obj.rebSim(is, ip).aum_wTC(end,1);
                end
                
            end % end loop simulations
        end % simulationRebalance
        
        % Extract the simulation statistics
        function [simStatsTab, TRC, CEC, TOT] = extractStatistics(obj)
            
            TRC = mean([obj.rebStats.TRC],1);
            CEC = mean([obj.rebStats.CEC],1);
            TOT = mean([obj.rebStats.TOT],1);
            
            simStatsTab = table(obj.strNames', TRC', CEC', TOT', ...
                'VariableNames', {'Reb Strategy', 'TC', 'CEC', 'TOT'});
           
        end

    end % methods public
    
    methods (Access = protected)
        function propgrp = getPropertyGroups(obj)
            if ~isscalar(obj)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
            else
                if  ~isempty(obj.QAprx) && (sum(sum(abs(tril(obj.QAprx,-1)))) < 100*eps && sum(sum(abs(triu(obj.QAprx,1)))) < 100*eps)
                    diagOnly = 1;
                    diagParam = obj.QAprx(1,1);
                elseif ~isempty(obj.QAprx) && (sum(sum(tril(obj.QAprx,-1))) ~= 0 || sum(sum(triu(obj.QAprx,1))) ~= 0)
                    diagOnly = 2;
                else % QAprx has not been set
                    diagOnly = 0;
                end
                if diagOnly == 1 
                    propList = struct('funcName', obj.funcName, 'assetNames', obj.assetNames, 'dateCreated', obj.dateCreated, ...
                                      'optWts', obj.optWts, 'QAprx_OnlyDiag', diagOnly == 1, 'diagParam', diagParam, ...
                                      'avgTC', obj.avgTC, 'avgCEC', obj.avgCEC, 'numSims',obj.numSims, ...
                                      'numPeriods', obj.numPeriods, 'expRts', obj.expRts,   ...
                                      'tnsCst', obj.tnsCst, 'covMtx', obj.covMtx,  'lambdaAbs', obj.lambdaAbs,   ...
                                      'lambdaRel', obj.lambdaRel );
                elseif diagOnly == 2 
                    propList = struct('funcName', obj.funcName, 'assetNames', obj.assetNames, 'dateCreated', obj.dateCreated, ...
                                      'optWts', obj.optWts, 'QAprx_OnlyDiag', diagOnly == 2, 'QAprx', obj.QAprx, ...
                                      'avgTC', obj.avgTC, 'avgCEC', obj.avgCEC, 'numSims',obj.numSims, ...
                                      'numPeriods', obj.numPeriods, 'expRts', obj.expRts,  ...
                                      'tnsCst', obj.tnsCst, 'covMtx', obj.covMtx,  'lambdaAbs', obj.lambdaAbs,   ...
                                      'lambdaRel',obj.lambdaRel );
                else
                    propList = struct('funcName', obj.funcName, 'assetNames', obj.assetNames, 'dateCreated', obj.dateCreated, ...
                                      'curWts', obj.curWts, 'newWts', obj.newWts, 'optWts', obj.optWts, 'bentWts', obj.benWts, ...
                                      'tnsCst', obj.tnsCst, 'expRts', obj.expRts,  ...
                                      'covMtx', obj.covMtx,  'lambdaAbs', obj.lambdaAbs,   ...
                                      'lambdaRel',obj.lambdaRel );
                end
                propgrp = matlab.mixin.util.PropertyGroup(propList);
            end
        
            
        end % getProperties
    end % methods protected
    
end

% Mean-variance-tracking
function f = cfMVT(w, wb,  mu, lambda_abs, lambda_rel, cov_matrix)

    f = -mu*w + lambda_abs*w'*cov_matrix*w + lambda_rel*(w - wb)' * cov_matrix * (w - wb);

end