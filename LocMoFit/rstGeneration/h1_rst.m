function h1_rst(str)
    global fid
    str_len = length(str);
    str2 = repelem('=', str_len);
    str = [str newline str2];
    fprintf(fid, [str newline]);
end
