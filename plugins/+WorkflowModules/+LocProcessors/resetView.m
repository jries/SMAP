classdef resetView<interfaces.WorkflowModule;
    properties

        
    end
    methods
       function obj=resetView(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)
           
        end
        function output=run(obj,data,p)
             si=obj.getPar('sr_sizeRecPix');
              if ~isempty(obj.locData.loc)
                  if p.setview
                      maxx=p.xrange(2);minx=p.xrange(1);
                      maxy=p.yrange(2);miny=p.yrange(1);
                  else
                        mx=myquantilefast(obj.locData.loc.xnm,[0.9995,0.0005],100000);
                        maxx=mx(1);minx=mx(2);
                        my=myquantilefast(obj.locData.loc.ynm,[0.9995,0.0005],100000);
                        maxy=my(1);miny=my(2);
                  end
              else
                  disp('cannot find size of image, no reset')
                  return
              end
              obj.setPar('sr_pos',[(maxx+minx)/2 (maxy+miny)/2]);
              pixrec=round(max((maxx-minx)/si(1),(maxy-miny)/si(2)));
              obj.setPar('sr_pixrec',pixrec);
              output=data;
                    
        end
    end
end


function pard=guidef

pard.setview.object=struct('Style','checkbox','String','Set view to [min max] (nm)','Value',0);
pard.setview.position=[1,1];
pard.setview.Width=1;

pard.tx.object=struct('Style','text','String','X');
pard.tx.position=[1,2];
pard.tx.Width=.3;

pard.xrange.object=struct('Style','edit','String','0 55000');
pard.xrange.position=[1,2.3];
pard.xrange.Width=.7;

pard.ty.object=struct('Style','text','String','Y');
pard.ty.position=[1,3];
pard.ty.Width=.3;

pard.yrange.object=struct('Style','edit','String','0 55000');
pard.yrange.position=[1,3.3];
pard.yrange.Width=.7;


pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='resets view';
end