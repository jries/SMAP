classdef PushToFiji<interfaces.DialogProcessor
    % PushToFiji opens the reconstructed superresolution image in Fiji
    methods
        function obj=PushToFiji(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.showresults=false;
        end
        
        function out=run(obj,p)
            out=[];
            mij=openfiji(obj);
             title=obj.getPar('layer1_').ch_filelist.selection;
             outimage=makeoutputtif(obj,p);
%             switch p.outputformat.selection
%                 case 'rendered image with scalebar'
%                     srimage=obj.getPar('sr_image');
%                     outimage=uint8(srimage.image*255);
% %                        mij.createColor(title,outimage,true);
%                 case 'rendered image'
%                     srimage=obj.getPar('sr_image');
%                     outimage=uint8(srimage.composite*255);
% %                        mij.createColor(title,outimage,true);
%                 case 'layer 1-3 as RGB'
%                     s1=obj.locData.layer(1).images.srimage.image;
%                     sizes=size(s1);
%                     if length(sizes)>2
%                         s1=sum(s1,3);
%                     end
%                     outimage=zeros(sizes(1),sizes(2),3);
%                     outimage(:,:,1)=s1;
%                     if length(obj.locData.layer)>1
%                          s2=obj.locData.layer(2).images.srimage.image;
%                           if length(size(s2))>2
%                             s2=sum(s2,3);
%                           end
%                           if size(s2)==size(s1)
%                               outimage(:,:,2)=s2;
%                           end
%                     end
%                     if length(obj.locData.layer)>2
%                          s3=obj.locData.layer(3).images.srimage.image;
%                           if length(size(s3))>2
%                             s3=sum(s3,3);
%                           end
%                           if size(s3)==size(s1)
%                               outimage(:,:,3)=s3;
%                           end
%                     end
%                     outimage=uint8(outimage/max(outimage(:))*255);
% %                        mij.createColor(title,outimage,true);
%                         
%                 case 'layer1 as grayscale'
%                     s1=obj.locData.layer(1).images.srimage.image;
%                     sizes=size(s1);
%                     if length(sizes)>2
%                         s1=sum(s1,3);
%                     end
%                     outimage=s1;
%                     outimage=uint16(outimage/max(outimage(:))*(2^16));
% %                        mij.createImage(title,outimage,true);
%             end      
            img=copytoImagePlus(outimage);
            img.show;
        end
        function exitfiji(obj,a,b)
            mij=obj.getPar('IJM');
            if ~isempty(mij)
                obj.setPar('IJM',[]);
                mij.exit();
            end
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            obj.guihandles.outputformat.String=makeoutputtif;
        end
    end
end

function pard=guidef(obj)
pard.t1.object=struct('String','select output figure','Style','text');
pard.t1.position=[1,1];
pard.outputformat.object=struct('String',{{'rendered image with scalebar','rendered image','layer 1-3 as RGB','layer1 as grayscale'}},'Style','popupmenu');
pard.outputformat.position=[2,1];
pard.outputformat.Width=2;

           pard.addscalebar.object=struct('String','add scale bar','Style','checkbox');
            pard.addscalebar.position=[2,3];
            pard.addscalebar.Width=2;
% pard.exitfiji.object=struct('String','exit Fiji','Style','pushbutton','Callback',@obj.exitfiji);
% pard.exitfiji.position=[4,1];
% pard.exitfiji.Width=1;

pard.plugininfo.name='Open in Fiji';
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='opens the reconstructed superresolution image in Fiji';

end

