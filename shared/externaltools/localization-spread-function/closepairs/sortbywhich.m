function ind = sortbywhich(x,y,rmax,t,taumin,taumax)

arguments
    x       (1,:)       double
    y       (1,:)       double
    rmax    (1,1)       double
    t       (1,:)       double = []
    taumin  (1,1)       double = 0
    taumax  (1,1)       double = 0
end

if isempty(t) % only x and y matter
    % pick whichever has the wider range
    ind = double(max(x) - min(x)) < (max(y) - min(y)) + 1;
else 
    % get naive expected number of comparisons for each direction
    xprs = 2*rmax/(max(x) - min(x));
    yprs = 2*rmax/(max(y) - min(y));
    trange = max(t) - min(t);
    tprs = (taumax -  taumin)/(trange - taumin);

    prs = [xprs, yprs, tprs];
    ind = find(prs == min(prs), 1);

end

end
