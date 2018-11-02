function [x,lb,ub,msg] = checkbounds(xin,lbin,ubin,nvars)
%CHECKBOUNDS Verify that the bounds are valid with respect to initial point.
%
% This is a helper function.

%   [X,LB,UB,X,FLAG] = CHECKBOUNDS(X0,LB,UB,nvars) 
%   checks that the upper and lower
%   bounds are valid (LB <= UB) and the same length as X (pad with -inf/inf
%   if necessary); warn if too long.  Also make LB and UB vectors if not 
%   already. Finally, inf in LB or -inf in UB throws an error.

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.4.2 $  $Date: 2006/12/15 19:29:30 $

msg = [];
% Turn into column vectors
lb = lbin(:); 
ub = ubin(:); 
xin = xin(:);

lenlb = length(lb);
lenub = length(ub);
lenx = length(xin);

% Check maximum length
if lenlb > nvars
   warning('optimlib:checkbounds:IgnoringExtraLbs', ...
           'Length of lower bounds is > length(x); ignoring extra bounds.');
   lb = lb(1:nvars);   
   lenlb = nvars;
elseif lenlb < nvars
   lb = [lb; -inf*ones(nvars-lenlb,1)];
   lenlb = nvars;
end

if lenub > nvars
   warning('optimlib:checkbounds:IgnoringExtraUbs', ...
           'Length of upper bounds is > length(x); ignoring extra bounds.');
   ub = ub(1:nvars);
   lenub = nvars;
elseif lenub < nvars
   ub = [ub; inf*ones(nvars-lenub,1)];
   lenub = nvars;
end

% Check feasibility of bounds
len = min(lenlb,lenub);
if any( lb( (1:len)' ) > ub( (1:len)' ) )
   count = full(sum(lb>ub));
   if count == 1
      msg=sprintf(['Exiting due to infeasibility:  %i lower bound exceeds the' ...
            ' corresponding upper bound.'],count);
   else
      msg=sprintf(['Exiting due to infeasibility:  %i lower bounds exceed the' ...
            ' corresponding upper bounds.'],count);
   end 
end
% check if -inf in ub or inf in lb   
if any(eq(ub, -inf)) 
   error('optimlib:checkbounds:MinusInfUb', ...
         '-Inf detected in upper bound: upper bounds must be > -Inf.');
elseif any(eq(lb,inf))
   error('optimlib:checkbounds:PlusInfLb', ...
         '+Inf detected in lower bound: lower bounds must be < Inf.');
end

x = xin;
