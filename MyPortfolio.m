classdef MyPortfolio
    properties
        assetNames
        dateCreated
        curWts
        benWts
        tnsCst
        expRts
        covMtx
        prtRet
        prtRsk
    end
    methods (Access = public)
        
        % Constructor
        function obj = MyPortfolio()
            obj.dateCreated = date();
        end % constructor
        
        % Set properties
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
        
        % Calculate expected portfolio returns
        function obj = portfolioReturn(obj)
            
            if ~isempty(obj.curWts) && (length(obj.curWts) == length(obj.expRts))
                obj.prtRet = obj.curWts * obj.expRts';
            end
            
        end % portfolioReturn
        
        % Calculate expected portfolio risk
        function obj = portfolioRisk(obj)
            
            if ~isempty(obj.curWts) && (length(obj.curWts) == size(obj.covMtx,1))
                obj.prtRsk = sqrt(obj.curWts * obj.covMtx * obj.curWts');
            end
            
        end
    end  % methods - public
    
end % classdef Portfolio
            