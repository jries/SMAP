function imout=renderallSMAP(locs,obj,varargin)
%render from any function using SMAP settings
pin=parseinput(varargin);
%later extend to input LocalizationData
%now structure
plotfields=pin.plotfields;
if isempty(pin.plotfields)
    if isfield(locs,'x')
        plotfields={'x','y'};
    elseif isfield(locs,'xnm')
        plotfields={'xnm','ynm'};
    end
end
locs.x=locs.(plotfields{1});
locs.y=locs.(plotfields{2});

fields=union(union(renderSMAP, drawerSMAP),obj.inputParameters);
obj.inputParameters=fields;
p=obj.getAllParameters;
pl=obj.getLayerParameters;
for k=length(pl):-1:1
     ph=pl{1};
    if ~ph.layercheck
        continue
    end
   
    ph.sr_axes=[];
    if isempty(pin.position)
        position=[(min(locs.x)+max(locs.x))/2,(min(locs.y)+max(locs.y))/2,max(locs.x)-min(locs.x),max(locs.y)-min(locs.y)];
        %x,y,w,h
    else
        position=pin.position;
    end
        ph.rangex=[position(1) position(1)+position(3)];
        ph.rangey=[position(2) position(2)+position(4)];
    if ~isempty(pin.pixrec)
        ph.sr_pixrec=pin.pixrec;
    end
    images.srimage=renderSMAP(locs,ph,1);
    images.finalImages=drawerSMAP(images.srimage,ph);
    layer(k).images=images;
end
 imout=displayerSMAP(layer,ph);

end

function pv=parseinput(in)
p=inputParser;

p.addParameter('position',[],@isnumeric);
p.addParameter('indin',[],@isnumeric);
p.addParameter('plotfields',[],@isnumeric);
p.addParameter('pixrec',[],@isnumeric);

parse(p,in{:});
pv=p.Results;
% pv
end