function resize(ga, p)
    % ga: graphics array
    % p: portion
    for k = 1:length(ga)
        h = ga(k);
        h.Position(1:2) = h.Position(1:2)+h.Position(3:4).*(1-p)/2;
        h.Position(3:4) = h.Position(3:4).*p;
    end
end