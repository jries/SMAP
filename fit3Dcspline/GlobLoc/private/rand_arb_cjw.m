function out = rand_arb_cjw(dimensions, distribution, plotting)

    if nargin < 3
        plotting = 0;
    end % if
% 
%     load([distribution '.mat'])

    if size(dimensions,2) == 1
        dimensions = [dimensions, 1];
    end % if

    x = distribution(:,1);
    y = distribution(:,2);
    y = y./max(y);

    out = NaN(dimensions);

    for m = 1:size(out,1)
        for n = 1:size(out,2)
            while isnan(out(m,n))
                x_rand = rand(1)*(max(x)-min(x)) + min(x);
                y_rand = rand(1);
                y_funct = y(find(x./x_rand>=1,1));
                if y_rand < y_funct
                    out(m,n) = x_rand;
                end % if
            end % while
        end % for n
    end % for m

    if plotting == 1
        [y_out, x_out] = hist(out, 100);
        figure
            plot(x,y, '--')
        hold on
        plot(x_out, y_out./max(y_out))
    end % if

end % function
