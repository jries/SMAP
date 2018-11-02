classdef Loc2pos<interfaces.WorkflowModule
    properties
        filestruc;
        locs
    end
    methods
        function obj=Loc2pos(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,data,p)
            global intLoc2pos_ind2 intLoc2pos_locframes;
            intLoc2pos_ind2=1;
            obj.locs=obj.locData.getloc({'frame','xnm','ynm'},'position','roi','layer',find(obj.getPar('sr_layerson')));
          
            intLoc2pos_locframes=obj.locs.frame;
            obj.filestruc=obj.locData.files.file(1);

        end
        function datout=run(obj,data,p)
            global intLoc2pos_ind2 intLoc2pos_locframes

            
            lf=length(obj.locs.xnm);
            frame=data.frame;
                %find indices for same frame
                ind1=intLoc2pos_ind2;
                while ind1>0&&intLoc2pos_locframes(ind1)<frame && ind1<lf
                    ind1=ind1+1;
                end
                ind1=min(ind1,lf);
                intLoc2pos_ind2=ind1;
                if ~intLoc2pos_locframes(intLoc2pos_ind2)==frame %no localizatiaon in frame
                      datout=data;%.copy;
                     datout.data.x=[];%.set(maxout);
                    return
                end
                while intLoc2pos_ind2<=lf&&intLoc2pos_locframes(intLoc2pos_ind2)==frame
                    intLoc2pos_ind2=intLoc2pos_ind2+1;
                end
                intLoc2pos_ind2=intLoc2pos_ind2-1;
 
              [~,maxout]=nm2pixLoc(obj.locs.xnm(ind1:intLoc2pos_ind2),obj.locs.ynm(ind1:intLoc2pos_ind2),obj.filestruc.info.cam_pixelsize_um*1000,obj.filestruc.info.roi);
               maxout.frame=frame+0*maxout.y;
               datout=data;%.copy;
               datout.data=maxout;%.set(maxout);
%                obj.output(datout); 
                if intLoc2pos_ind2==lf
                    datout.eof=true;
                end
        end
    end
end




function pard=guidef

pard.plugininfo.type='WorkflowModule'; 
end
