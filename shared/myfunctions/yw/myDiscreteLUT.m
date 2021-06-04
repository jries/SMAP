function output = myDiscreteLUT(idx)
% Discrete used in Robin's paper on nature methods
    lut = {'#f80503', '#0606c3', '#fa8256', '#aaaaaa', '#05fa9d', '#8a2be0' '#050505' '#0ccacd', '#e7e444'};
    if length(idx) == 1
        output = lut{idx};
    else
        output = lut(idx);
    end
end
