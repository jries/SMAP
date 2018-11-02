function ig=anyfilter(v,mode,varargin)
ig=true(size(v));
switch mode
    case 'minmax'
        if length(varargin)==1
            if iscell(varargin{1})
                minx=varargin{1}{1};
                maxx=varargin{1}{2};
            else
                minx=varargin{1}(1);
                maxx=varargin{1}(2);
            end
        else
            minx=varargin{1};
            maxx=varargin{2};
        end
        if ~isempty(v)
            ig=v>=minx&v<=maxx;    
        else
            ig=[];
        end
    case 'inlist'
        if ~isempty(v)
            ig=false(size(v));
            for k=1:length(varargin{1})
                ig=ig|v==varargin{1}(k);
            end 
        end
        
end