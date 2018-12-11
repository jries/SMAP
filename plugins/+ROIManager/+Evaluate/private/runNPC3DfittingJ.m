function out = runNPC3DfittingJ(obj, p)
    % Set bounds of parameters
    locs = obj.getLocs({'xnmrot','ynmrot', 'znm'},'size',obj.P.par.se_siteroi.content/2,'grouping','ungrouped', 'layer', 1);

    % Build the model/template
    img3d = obj.getPar('img3d');
    xcor3d =  obj.getPar('xcor3d');
    ycor3d =  obj.getPar('ycor3d');
    zcor3d =  obj.getPar('zcor3d');
    if isempty(img3d)
        img3d = generate3dRing([]);
        [xcor3d, ycor3d, zcor3d]=meshgrid(1:length(img3d), 1:length(img3d), 1:length(img3d));
        obj.setPar('img3d', img3d);
        obj.setPar('xcor3d', xcor3d);
        obj.setPar('ycor3d', ycor3d);
        obj.setPar('zcor3d', zcor3d);
    end
  p.sumOffset=1;
    % Objective function
    objfun =@(x) -sum(log10(NPC3dJ(locs.xnmrot, locs.ynmrot, locs.znm, x(1), x(2), x(3), x(4), x(5), x(6) ,'gridX', xcor3d, 'gridY',ycor3d, 'gridZ',zcor3d, 'img3d',img3d, 'roisize', 500, 'SumOffset', p.sumOffset)));
%      objfun =@(x) -sum((NPC3d(locs.xnmrot, locs.ynmrot, locs.znm, x(1), x(2), x(3), x(4), x(5), x(6) ,'gridX', xcor3d, 'gridY',ycor3d, 'gridZ',zcor3d, 'img3d',img3d, 'roisize', 500, 'SumOffset', p.sumOffset)));
    options=optimset('TolX',1e-7);
    dstart=60;
    st=[double(median(locs.xnmrot)) double(mean(locs.xnmrot)) double(mean(locs.znm))-dstart dstart 0 0 ];
    x=st;
    
%           NPC3d(locs.xnmrot, locs.ynmrot, locs.znm, x(1), x(2), x(3), x(4), x(5), x(6), 'gridX', xcor3d, 'gridY',ycor3d, 'gridZ',zcor3d, 'img3d',img3d, 'roisize', 500, 'SumOffset', p.sumOffset, 'plot',1)
          
   fittedVals= fminsearch(objfun,st,options);

%     if 0
%         errors = 20 * rand([30 1]) - 10;
%         intDist = 40 + errors; 
%         boundRange = p.ub - p.lb;
%         boundMean = (p.ub + p.lb)/2;
%         boundRange = boundRange*0.5;
%         intRand = rand([30 6]);
%         intRand = intRand .* boundRange;
%         intRand = intRand + boundMean;
%         intRand = [intRand(:,1:3) intDist intRand(:,5:end)];
%         options = optimoptions('particleswarm', 'InitialSwarmMatrix', intRand);
%     end
%         options = optimoptions('particleswarm');
%         fittedVals = particleswarm(objfun, 6, p.lb, p.ub, options);
        out.fitted = fittedVals;
        x = fittedVals;
        [~,~,~,~,~,out.locxyz] = NPC3dJ(locs.xnmrot, locs.ynmrot, locs.znm, x(1), x(2), x(3), x(4), x(5), x(6), 'gridX', xcor3d, 'gridY',ycor3d, 'gridZ',zcor3d, 'img3d',img3d, 'roisize', 500, 'SumOffset', p.sumOffset, 'plot',[]);
    if obj.display
        ax=obj.setoutput('fit');
        NPC3dJ(locs.xnmrot, locs.ynmrot, locs.znm, x(1), x(2), x(3), x(4), x(5), x(6), 'gridX', xcor3d, 'gridY',ycor3d, 'gridZ',zcor3d, 'img3d',img3d, 'roisize', 500, 'SumOffset', p.sumOffset, 'plot',ax)
    end
end