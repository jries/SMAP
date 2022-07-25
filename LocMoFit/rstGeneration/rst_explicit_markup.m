function rst_explicit_markup(markup, ref, varargin)
    global fid
    str = ['.. ' markup ':: ' ref];
    for k = 1:length(varargin)/2
        par = varargin{(k-1)*2+1};
        val = varargin{k*2};
        str2 = [newline char(9) ':' par ':' val];
        str = [str str2];
    end
    fprintf(fid, [str newline]);
end