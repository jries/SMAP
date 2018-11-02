function varargout = xerrorbar(varargin)
%XERRORBAR Error bar plot for x - values.
%   XERRORBAR(X,Y,L,U) plots the graph of vector X vs. vector Y with
%   error bars specified by the vectors L and U.  L and U contain the
%   lower and upper error ranges for each point in X.  Each error bar
%   is L(i) + U(i) long and is drawn a distance of U(i) left of and L(i)
%   right of the points in (X,Y).  The vectors X,Y,L and U must all be
%   the same length.  If X,Y,L and U are matrices then each column
%   produces a separate line.
%
%   XERRORBAR(X,Y,E) plots X with error bars [X-E X+E]
%
%   XERRORBAR(AX,...) uses the specified axis, instead of gca.
%
%   XERRORBAR(...,'LineSpec') uses the color and linestyle specified by
%   the string 'LineSpec'.  It applies these specifications to all the
%   lineseries specified and their error bars.
%
%   [HH] = XERRORBAR(...) returns a vector of handles to lineseries
%   objects.
%
%    For Example, we might have a dataseries with constant measurement
%    errors in x.
%       x = 0:10;
%       y = sin(x*pi/5);
%       e = .1*ones(size(x));
%       xerrorbar(x,y,e)
%     
%     
%   ---
%   MFILE:   xerrorbar.m
%   VERSION: 1.0 (2012/01/12) 
%   MATLAB:  7.8.0.347 (R2009a)
%   AUTHOR:  Anthony Bathgate
%   CONTACT: tony.bathgate@usask.ca
%
%   ADDITIONAL NOTES:
%     - I think it works just like errorbars, except it applies markers to 
%     the ends of the error bars as well as the data series itself.
%     
%
%   REVISIONS:
%   1.0      Released. (2011/09/15)
%   
%   DISCLAIMER:
%   xerrorbar.m is provided "as is" without warranty of any kind, 
%   under the revised BSD license.
%
%   Copyright (C) 2012 Tony Bathgate, University of Saskatchewan, ISAS.


    % Check I/O
    error(nargchk(3,inf,nargin,'struct'));
    error(nargoutchk(0,1,nargout,'struct'));
    
    % Figure out axes
    if( isa(handle(varargin{1}),'hg.axes') )
        ax      = varargin{1};
        args    = varargin(2:end);
        nargs   = nargin-1;
    else
        ax      = gca;
        args    = varargin;
        nargs   = nargin;
    end
    holdoff = ~ishold(ax);
    
    % Check vectors, and arrange them properly
    if(~samedims(args{1},args{2}) || ~samedims(args{2},args{3}))
        error('Argument dimension mismatch');
    end
    xd      = size(args{1});
    if( xd(1)==1 )
        args{1} = args{1}';
        args{2} = args{2}';
        args{3} = args{3}';
    end
    X       = args{1};
    Y       = args{2};
    L       = args{3};
    
    % If they provided a U~=L then use it. Otherwise use L.
    if( nargs>3 )
        if( isnumeric(args{4}) )
            hasu = true;
            if(xd(1)==1)
                args{4} = args{4}';
            end
            if( ~samedims(args{3},args{4}) )
                error('Argument dimension mismatch');
            end
            U = args{4};
        else
            hasu = false;
            U = L;
        end
    else
        hasu = false;
        U = L;
    end
    
    % Assemble the Line specifications
    lcol    = Lines;
    xd      = size(args{1});
    lspec   = cell(xd(2),1);
    for sctr = 1:xd(2)
        switch( nargs )
            case 3
                % Give each line a different color by default
                lspec{sctr} = {'Color',lcol(mod(sctr,length(lcol)),:)};
            case 4
                if( hasu )
                    % Give each line a different color by default
                    lspec{sctr} = {'Color',lcol(mod(sctr,length(lcol)),:)};
                else
                    % They've given a LineSpec so use it
                    lspec{sctr} = args(4);
                end
            otherwise
                % They've given more LineSpecs, so use them
                if( hasu )
                    lspec{sctr} = args(5:end);
                else
                    lspec{sctr} = args(4:end);
                end
        end
    end
    
    % Plot each line series
    hh = NaN(xd(2),1);
    hh(1) = xeb1(ax,X(:,1),Y(:,1),...
        L(:,1),U(:,1),lspec{1});
    if(holdoff)
        hold on;
    end
    for sctr = 2:xd(2)
        hh(sctr) = xeb1(ax,X(:,sctr),Y(:,sctr),...
            L(:,sctr),U(:,sctr),lspec{sctr});
    end

    % They want the handles
    if( nargout==1 )
        varargout{1} = hh;
    end

    % Reset the plot
    if(holdoff)
       hold off; 
    end
    
end

% Construct a single series of points that represents the entire lineseries
% (including error bars);
function [hh] = xeb1(ax,X,Y,L,U,lspec)
    ebtw = (max(Y)-min(Y))/100;
    
    xlong = X;
    ylong = Y;
    for pctr = 1:length(X)
        xlong = [xlong;NaN;X(pctr)-L(pctr);X(pctr)+U(pctr);];
        ylong = [ylong;NaN;Y(pctr);Y(pctr);];
        xlong = [xlong;NaN;X(pctr)-L(pctr);X(pctr)-L(pctr);];
        ylong = [ylong;NaN;Y(pctr)-ebtw;Y(pctr)+ebtw;];
        xlong = [xlong;NaN;X(pctr)+U(pctr);X(pctr)+U(pctr);];
        ylong = [ylong;NaN;Y(pctr)-ebtw;Y(pctr)+ebtw;];
    end
    
    hh = plot(ax,xlong,ylong,lspec{1:end});
end

% Check that all the dimensions match
function v = samedims(A,B)
    sa = size(A);
    sb = size(B);
    v = true;
    v = v && length(sa)==length(sb);
    if(v)
        for c = 1:length(sa)
            v = v&&sa(c)==sb(c);
        end
    end
end
