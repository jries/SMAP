function [cutoff, bin_edges, bin_n, curve] = getOneRingCutoff(ringDist)
    % Determine the cutoff between one-ring and dual-ring NPCs.
    
    bin_edges = floor(min(ringDist)):5:ceil(max(ringDist));
    bin_n = histcounts(ringDist,bin_edges);
    bin_center = movmean(bin_edges,2,'Endpoints', 'discard');

    % fit two gauss
    curve = fit(bin_center', bin_n','gauss2','StartPoint', [15,max(bin_n)/2,5,80,max(bin_n),10],'Lower', [0, 10, 0, 0,20,0],'Upper', [Inf, 30, 30, Inf,80,30],'Robust', 'LAR');
    gauss = @(x,a,b,c) a.*exp(-((x-b)./c).^2);

    % get the intersection between the two gauss
    cutoff = fzero(@(x) gauss(x,curve.a1,curve.b1,curve.c1) - gauss(x,curve.a2,curve.b2,curve.c2), (curve.b1+curve.b2)/2);
end