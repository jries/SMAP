classdef selectCoordinates2Color<interfaces.SEEvaluationProcessor
    properties
        hax1
        hax2
        hax3
        line
    end
    methods
        function obj=selectCoordinates2Color(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
            obj.line=p.lineselect.selection;
            switch obj.line
                case 'line1'
                    cellline='line2';
                otherwise
                    cellline='line1';
            end
            %  use shift_xy to fit with offset
            ax1=obj.setoutput('positions',true);
            hp=ax1.Parent;
%             delete(hp.Children);
%             ax1=axes(hp);
            ax2=axes(hp);
            ax3=axes(hp);
            ax4=axes(hp);
            ax5=axes(hp);
            ax6=axes(hp);
            subplot(2,3,1,ax1);
            subplot(2,3,2,ax2);
            subplot(2,3,3,ax3);
            subplot(2,3,4,ax4);
            subplot(2,3,5,ax5);
            subplot(2,3,6,ax6);
            
            site=obj.site;
            
%             cellind=obj.locData.SE.indexFromID(obj.locData.SE.cells,obj.site.info.cell);
            cell=obj.locData.SE.getcell(obj.site.info.cell);
            im1=site.image.layers(1).images.finalImages;
            if length(site.image.layers)<2
                ind=1;
            else
                ind=2;
            end
            im2=site.image.layers(ind).images.finalImages;
            
            range=[-p.se_sitefov p.se_sitefov]/2;
            imagesc(ax1,im1.rangex/1000,im1.rangey/1000,im1.image);
            axis(ax1,'equal');
            imagesc(ax2,im2.rangex/1000,im2.rangey/1000,im2.image);
            axis(ax2,'equal');
            imagesc(ax3,site.image.rangex,site.image.rangey,site.image.image);
            axis(ax3,'equal');
            pl=site.annotation.(obj.line).pos;
            if sum(pl)==0
                xp=mean(site.image.rangex);
                yp=mean(site.image.rangey);
                pl=[xp, yp;xp,yp];
                obj.site.annotation.(obj.line).pos=pl;
            end
            obj.hax1=impoint(ax1,pl(1,:));
                
            
            %direction of cell
            pcell=site.annotation.(cellline).pos;
            vcell=(pcell(2,:)-pcell(1,:));vcell=vcell/norm(vcell);
            
            cim1=cell.image.layers(1).images.finalImages;
            cim2=cell.image.layers(ind).images.finalImages;
%             tA=eye(3);
            dx=vcell(1)*p.shiftval/p.se_cellpixelsize;dy=vcell(2)*p.shiftval/p.se_cellpixelsize;
            cim2sp=imtranslate(cim2.image,[dx,dy]);
%              tA(3,1)=-vcell(1)*p.shiftval/p.se_cellpixelsize;tA(3,2)=-vcell(2)*p.shiftval/p.se_cellpixelsize;
            cim2sm=imtranslate(cim2.image,[-dx,-dy]);
            imagesc(ax4,cell.image.rangex,cell.image.rangey,cell.image.image)
            axis(ax4,'equal');
            imagesc(ax5,cell.image.rangex,cell.image.rangey,cim1.image+cim2sp)
            axis(ax5,'equal');
            imagesc(ax6,cell.image.rangex,cell.image.rangey,cim1.image+cim2sm)
            axis(ax6,'equal')
     
            
            obj.hax1.addNewPositionCallback(@(p) obj.posconstraint(p,1));
            obj.hax2=impoint(ax2,pl(2,:));
            obj.hax2.addNewPositionCallback(@(p) obj.posconstraint(p,2));
            obj.hax3=imline(ax3,pl);
            obj.hax3.addNewPositionCallback(@(p) obj.posconstraint(p,3));
            
            

            pcell(1,:)=pcell(1,:)-p.se_cellfov*vcell/1000;
            pcell(2,:)=pcell(2,:)+p.se_cellfov*vcell/1000;
            hcl=imline(ax4,pcell);
            
            
            out=[];

        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function posconstraint(obj,pos,caller)
            obj.line
            linepos=obj.site.annotation.(obj.line).pos;
            
            switch caller
                case 1
                    linepos(1,:)=pos;
                case 2
                    linepos(2,:)=pos;
                case 3
                    linepos=pos;
            end
               
            obj.hax1.setPosition(linepos(1,:));
            obj.hax2.setPosition(linepos(2,:));
            obj.hax3.setPosition(linepos);
            obj.site.annotation.(obj.line).pos=linepos;
            
%             obj.site.annotation.line1.value=sqrt(sum((linepos(2,:)-linepos(1,:)).^2));
%             obj.site.annotation.line1.length=obj.site.annotation.line1.value;
        end
    end
end

function pard=guidef
pard.inputParameters={'se_cellpixelsize','numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer2_'};

pard.lineselect.object=struct('Style','popupmenu','String',{{'line1','line2'}});
pard.lineselect.position=[1,1];
pard.lineselect.Width=4;

pard.shiftvalt.object=struct('Style','text','String','mean shift nm');
pard.shiftvalt.position=[2,1];
pard.shiftvalt.Width=2;

pard.shiftval.object=struct('Style','edit','String','0');
pard.shiftval.position=[2,3];
pard.shiftval.Width=2;

pard.plugininfo.type='ROI_Evaluate';

end


