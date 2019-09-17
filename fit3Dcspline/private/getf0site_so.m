function [f0,PSFx0,PSFy0]=getf0site_so(locs,p)
    p.ploton=0;
    if isfield(locs,'znm')&&~isempty(locs.znm)
         zas=getf0Z(locs,p);
    else
    [zas,zn]=stackas2z(locs.PSFxpix,locs.PSFypix,locs.frame,locs.phot,p);
    end
    if isnan(zas)
        PSFx0=NaN;
        PSFy0=NaN;
    else
        ind=find(locs.frame<=zas,1,'last');
        if isempty(ind)
            ind=1;
        end
    PSFx0=locs.PSFxpix(ind);
    PSFy0=locs.PSFypix(ind);
    end
    % drawnow
    f0=zas;

end