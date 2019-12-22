classdef Costfunction < MyPortfolio & matlab.mixin.CustomDisplay

    properties
        funcName
        newWts
        optWts
        lambdaAbs
        lambdaRel
        funcValCur
        funcValNew
        funcValOpt
        paidTC
        exitFlag
    end
    methods (Access = public)
        % Constructor
        function obj = Costfunction()
            
            obj@MyPortfolio();
            
        end % constructor Costfunc
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
        
        % Optimize costfunctions
        function obj = optimizeCostfunction(obj, getOptimal)
            
            % Parameters
            numAssets = length(obj.curWts);
            
            % fmincon options
            options = optimoptions('fmincon');
            options = optimoptions(options,'Display', 'notify', ...
                        'ConstraintTolerance', 10*eps, ...
                        'OptimalityTolerance', 10*eps);
                    
            % Contraint parameters and initial guess
            dw0 = zeros(2*numAssets,1);
            A = [];
            b = [];
            Aeq = [ones(1,numAssets) -ones(1,numAssets)];
            beq = 0;
            lb = zeros(2*numAssets,1);
            ub = [1 - obj.curWts'; obj.curWts'];
            nonlincons = [];

            switch obj.funcName
                case 'MV'
                    if ~getOptimal % Optimize with transaction costs
                        [dw, ~, obj.exitFlag] = fmincon(@(dw)cfMVT_wTC(obj, dw, obj.curWts', ...
                                            zeros(numAssets,1), obj.expRts, ...
                                            obj.tnsCst, obj.lambdaAbs, 0, obj.covMtx), ...
                                    dw0, A, b, Aeq, beq, lb, ub, nonlincons, options);
                        obj.newWts = obj.curWts + (dw(1:end/2) - dw(end/2+1:end))';
                        obj.funcValCur = cfMVT(obj, obj.curWts', zeros(numAssets,1), obj.expRts, ...
                                             obj.lambdaAbs, 0, obj.covMtx);
                        obj.funcValNew = cfMVT(obj, obj.newWts', zeros(numAssets,1), obj.expRts, ...
                                             obj.lambdaAbs, 0, obj.covMtx);
                        obj.paidTC = obj.tnsCst * abs(obj.curWts - obj.newWts)';
                    else % Optimize without transaction costs -> OPTIMAL ALLOCATION
                        [dw, ~, obj.exitFlag] = fmincon(@(dw)cfMVT_wTC(obj, dw, obj.curWts',  ...
                                            zeros(numAssets,1), obj.expRts,  ...
                                    zeros(size(obj.tnsCst)), obj.lambdaAbs, 0, obj.covMtx), dw0, ...
                                     A, b, Aeq, beq, lb, ub, nonlincons, options);            
                        obj.optWts = obj.curWts + (dw(1:end/2) - dw(end/2+1:end))';
                        obj.funcValCur = cfMVT(obj, bj.curWts', zeros(numAssets,1), obj.expRts, ...
                                             obj.lambdaAbs, 0, obj.covMtx);
                        obj.funcValOpt = cfMVT(obj, obj.optWts', zeros(numAssets,1), obj.expRts, ...
                                             obj.lambdaAbs, 0, obj.covMtx);
                    end
                case 'MVT'
                    if ~getOptimal % Optimize with transaction costs
                        [dw, ~, obj.exitFlag] = fmincon(@(dw)cfMVT_wTC(obj, dw, obj.curWts', ...
                                            obj.benWts', obj.expRts, obj.tnsCst, obj.lambdaAbs, ...
                                    obj.lambdaRel, obj.covMtx), dw0, A, b, Aeq, beq, lb, ub, ...
                                    nonlincons, options);
                        obj.newWts = obj.curWts + (dw(1:end/2) - dw(end/2+1:end))';
                        obj.funcValCur = cfMVT(obj, obj.curWts', obj.benWts', obj.expRts,  ...
                                            obj.lambdaAbs, obj.lambdaRel, obj.covMtx);
                        obj.funcValNew = cfMVT(obj, obj.newWts', obj.benWts', obj.expRts,  ...
                                            obj.lambdaAbs, obj.lambdaRel, obj.covMtx);
                        obj.paidTC = obj.tnsCst * abs(obj.curWts - obj.newWts)';
                    else % Optimize without transaction costs
                        [dw, ~, obj.exitFlag] = fmincon(@(dw)cfMVT_wTC(obj, dw, obj.curWts', ...
                                            obj.benWts', obj.expRts, zeros(size(obj.tnsCst)),  ...
                                    obj.lambdaAbs, obj.lambdaRel, obj.covMtx), dw0, A, b, Aeq, ... 
                                    beq, lb, ub, nonlincons, options);
                        obj.optWts = obj.curWts + (dw(1:end/2) - dw(end/2+1:end))';
                        obj.funcValCur = cfMVT(obj, obj.curWts', obj.benWts', obj.expRts, ...
                                            obj.lambdaAbs, obj.lambdaRel, obj.covMtx);
                        obj.funcValOpt = cfMVT(obj, obj.optWts', obj.benWts', obj.expRts, ...
                                            obj.lambdaAbs, obj.lambdaRel, obj.covMtx);
                    end
                case 'QAprx' % We use MVT for simplicity MV = MVT with lambdaRel = 0
                    if ~getOptimal % Optimize with transaction costs and QAprx
                        [dw, ~, obj.exitFlag] = fmincon(@(dw)cfMVT_wTC_wQAprx(obj, dw, ...
                                           obj.curWts', obj.benWts', obj.optWts', obj.expRts, ... 
                                           obj.tnsCst, obj.lambdaAbs, obj.lambdaRel, obj.covMtx, ...
                                    obj.QAprx), dw0, A, b, Aeq, beq, lb, ub, nonlincons, options);
                            
                        obj.newWts = obj.curWts + (dw(1:end/2) - dw(end/2+1:end))';
                        obj.funcValCur = cfMVT(obj, obj.curWts', obj.benWts', obj.expRts, ...
                                            obj.lambdaAbs, obj.lambdaRel, obj.covMtx);
                        obj.funcValNew = cfMVT(obj, obj.newWts', obj.benWts', obj.expRts, ...
                                            obj.lambdaAbs, obj.lambdaRel, obj.covMtx);
                        obj.paidTC = obj.tnsCst * abs(obj.curWts - obj.newWts)';
                    else % Optimize without transaction costs and without QAprx
                        [dw, ~, obj.exitFlag] = fmincon(@(dw)cfMVT_wTC_wQAprx(obj, dw, ...
                                           obj.curWts', obj.benWts', obj.optWts', obj.expRts, ... 
                                           zeros(size(obj.tnsCst)), obj.lambdaAbs, obj.lambdaRel, ...
                                    obj.covMtx, zeros(size(obj.QAprx))), dw0, A, b, Aeq, beq, ...
                                    lb, ub, nonlincons, options);
                        obj.optWts = obj.curWts + (dw(1:end/2) - dw(end/2+1:end))';
                        obj.funcValCur = cfMVT(obj, obj.curWts', obj.benWts', obj.expRts, ...
                                            obj.lambdaAbs, obj.lambdaRel, obj.covMtx);
                        obj.funcValOpt = cfMVT(obj, obj.optWts', obj.benWts', obj.expRts, ...
                                            obj.lambdaAbs, obj.lambdaRel, obj.covMtx);
                    
                    end
            end
            
            % Exitflag check
            if ~(obj.exitFlag == 1 || obj.exitFlag == 2)
                fprintf("Failure to optimize -> 42 wasn't the answer apperently:%s\n", int2str(obj.exitFlag))
            end
            
        end % optimizeCostfunction
        
        % Mean-variance-tracking
        function f = cfMVT(obj, w, wb,  mu, lambda_abs, lambda_rel, cov_matrix)

            f = -mu*w + lambda_abs*w'*cov_matrix*w + lambda_rel*(w - wb)' * cov_matrix * (w - wb);

        end
        
        % Mean-variance-tracking with tranaction costs
        function f = cfMVT_wTC(obj, dw, wc, wb, mu, tnsCst, lambdaAbs,  lambdaRel, covMtx)

            % Split weights into positive and negative parts
            u = dw(1:end/2);
            v = dw(end/2+1:end);

            f = -mu*(wc +(u - v)) + lambdaAbs*(wc + (u - v))'*covMtx*(wc + (u - v)) + ...
                + lambdaRel*(wc - wb + (u - v))'*covMtx*(wc - wb +(u - v)) + tnsCst*(u + v);

        end % cfMVT_wTC

        % Mean-variance-tracking with tranaction costs and QAprx
        function f = cfMVT_wTC_wQAprx(obj, dw, wc, wb, wo, mu, tnsCst, lambdaAbs,  lambdaRel, covMtx, QAprx)

            % Split weights into positive and negative parts
            u = dw(1:end/2);
            v = dw(end/2+1:end);

            f = -mu*(wc +(u - v)) + lambdaAbs*(wc + (u - v))'*covMtx*(wc + (u - v)) + ...
                + lambdaRel*(wc - wb + (u - v))'*covMtx*(wc - wb +(u - v)) + tnsCst*(u + v) + ...
                ((wc + (u - v)) - wo)'*QAprx*((wc + (u - v)) - wo);

        end 
        
    end
    methods (Access = protected)
        function propgrp = getPropertyGroups(obj)
            
            if ~isscalar(obj)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
            else
            
                propList = struct('assetNames', obj.assetNames, 'dateCreated', obj.dateCreated, ...
                                  'curWts', obj.curWts, 'newWts', obj.newWts, 'optWts', obj.optWts, ...
                                  'benWts', obj.benWts, 'tnsCst', obj.tnsCst, 'expRts', obj.expRts, ...
                                  'covMtx', obj.covMtx, 'funcName', obj.funcName, 'lambdaAbs', obj.lambdaAbs,   ...
                                  'lambdaRel',obj.lambdaRel, 'funcValCur', obj.funcValCur,  ...
                                  'funcValNew',obj.funcValNew, 'funcValOpt',obj.funcValOpt, 'paidTC',obj.paidTC,  ...
                                  'exitFlag', obj.exitFlag);
                propgrp = matlab.mixin.util.PropertyGroup(propList);
            end
        
            
        end % getProperties
    end % methods protected
    
end % classdef Costfunc


