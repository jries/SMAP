function updateROI(g,varargin)
    % 'pos': position in um
    % 'rot': rotation in angle
    % 'width': in nm
    % 'length': in nm
    
    %% Get original parameters
    h = g.getPar('sr_roihandle');
    if isempty(h)
        warning('ROI does not exist. Please call createROI() first.')
        return
    end
    tips = getPosition(h);
    rot = lineAngle(tips);
    dist = tipDist(tips);
    
    % use the original parameters in default
    p = inputParser;
    p.addParameter('pos', mean(tips));
    p.addParameter('rot', rot);
    p.addParameter('width', g.getPar('linewidth_roi'));
    p.addParameter('length', dist*1000);
    parse(p,varargin{:})
    p = p.Results;
    
    %% Update ROI
    pos = p.pos;
    rot = p.rot;
    linewidth_roi = p.width;
    linelength_roi = p.length;
    prepos = [0 -linelength_roi/2;0 linelength_roi/2]./1000;
    [prepos(:,1),prepos(:,2)] = rotcoord(prepos(:,2),prepos(:,1),deg2rad(rot));
    pos = prepos+pos;

    setPosition(h,pos)
    g.setPar('linewidth_roi',linewidth_roi);
    notify(g.P,'sr_render')
end