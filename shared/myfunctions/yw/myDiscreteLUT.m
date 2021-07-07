function output = myDiscreteLUT(idx, varargin)
% Discrete used in Robin's paper on nature methods
    inp = inputParser;
    inp.addParameter('type', 'hex')
    inp.parse(varargin{:})
    inp = inp.Results;
    
    lut = {'#f80503', '#0606c3', '#fa8256', '#aaaaaa', '#05fa9d', '#8a2be0' '#050505' '#0ccacd', '#e7e444'};
    switch inp.type
        case 'hex'
            if length(idx) == 1
                output = lut{idx};
            else
                output = lut(idx);
            end
        case 'rgb'
            output = hex2rgb(lut(idx));
    end
    
end
