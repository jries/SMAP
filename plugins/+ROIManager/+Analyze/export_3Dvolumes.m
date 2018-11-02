classdef export_3Dvolumes<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=export_3Dvolumes(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        
        function out=run(obj,p)  
            sites=obj.SE.sites;
            if p.export_selected
                sev=obj.getPar('se_viewer');
                pv=sev.getAllParameters;
                selectedsites=pv.sitelist.Value;
            else
              
                selectedsites=1:length(sites);
            end
            mainfile=obj.getPar('mainfile');
            path=fileparts(mainfile);
            prefix='img_.em';
            [fileout,path]=uiputfile([path filesep prefix]);
            
            if ~fileout
                return
            end
            [~,f]=fileparts(fileout);
             layers=find(obj.getPar('sr_layerson'));
             pixelsize=p.pixrec;
             if length(pixelsize)<2
                pixelsize(2)=pixelsize(1);
            end
            if length(pixelsize)<3
                pixelsize(3)=pixelsize(1);
            end
            sizerec=p.sizerec;
            if length(sizerec)<2
                sizerec(2)=sizerec(1);
            end
            if length(sizerec)<3
                sizerec(3)=sizerec(1);
            end
            
             xrange=[-1 1]/2*sizerec(1)*pixelsize(1);
             yrange=[-1 1]/2*sizerec(2)*pixelsize(2);
             zrange=[-1 1]/2*sizerec(3)*pixelsize(3);
             
             positions=zeros(length(selectedsites),5);
            for k=selectedsites
                site=sites(k);
                for l=1:length(layers)
                    [locs,indloc]=obj.locData.getloc({'xnm','ynm','znm','numberInGroup'},'position',sites(k),'layer',layers(l));
                    zpos=median(locs.znm);
                    xh=locs.xnm-sites(k).pos(1);yh=locs.ynm-sites(k).pos(2);zh=locs.znm-zpos;
                    if p.blinkcoded
                        [imouth,dxo,dyo,dzo]=myhist3(xh,yh,zh,pixelsize,xrange,yrange,zrange,locs.numberInGroup);
                    else
                     imouth=renderhist3D(xh,yh,zh,xrange,yrange,zrange,pixelsize);
                    end
                    filen=[f strrep(site.name,'.','_') '_' int2str(l)];
                    switch p.format.selection
                    case '.em'

                        fhere= [filen '.em'];
                        tom_emwrite([path filesep fhere],imouth); 
                    otherwise
                        disp('format not implemented');
                    end
                end
                        
               
               positions(k,1:2)=sites(k).pos(1:2);
               positions(k,3)=zpos;
               positions(k,4)=sites(k).ID;
               positions(k,5)=k;
            end
            csvwrite([path filesep 'positions.txt'],positions);
            
            out=0;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function imouth=renderhist3D(x,y,z,xrange,yrange,zrange,pixelsize)
format='uint8';
if length(pixelsize)<2
    pixelsize(2)=pixelsize(1);
end
if length(pixelsize)<3
    pixelsize(3)=pixelsize(1);
end
 xh=ceil((x-xrange(1))/pixelsize(1));
 yh=ceil((y-yrange(1))/pixelsize(2));
 zh=ceil((z-zrange(1))/pixelsize(3));
 sizeV=ceil([(xrange(2)-xrange(1))/pixelsize(1) (yrange(2)-yrange(1))/pixelsize(2) (zrange(2)-zrange(1))/pixelsize(3)]);   
indg=xh>0&xh<=sizeV(1) & yh>0&yh<=sizeV(2) & zh>0&zh<=sizeV(3);
    
    linind=sub2ind(sizeV,xh(indg),yh(indg),zh(indg));
    hl=cast(histc(linind,1:max(linind)),format);
%     hl=cast(histcounts(linind,1:max(linind)),format);
%     disp(['maxcounts:' num2str(max(hl))])
    uind=unique(linind);
    imouth=zeros(sizeV,format);
    imouth(uind)=hl(hl>0);
end


function pard=guidef

pard.export_selected.object=struct('String','export only selected sites','Style','checkbox');
pard.export_selected.position=[1,1];
pard.export_selected.Width=2;

pard.format.object=struct('String',{{'.em','.tif'}},'Style','popupmenu');
pard.format.position=[2,1];
pard.format.Width=2;


pard.pixrect.object=struct('String','pixelsize (nm) [x y z]','Style','text');
pard.pixrect.position=[3,1];
pard.pixrect.Width=1;

pard.pixrec.object=struct('String','7','Style','edit');
pard.pixrec.position=[3,2];
pard.pixrec.Width=1;

pard.sizerect.object=struct('String','size (pixels) [x y z]','Style','text');
pard.sizerect.position=[4,1];
pard.sizerect.Width=1;

pard.sizerec.object=struct('String','50','Style','edit');
pard.sizerec.position=[4,2];
pard.sizerec.Width=1;

pard.blinkcoded.object=struct('String','blink coded for grouped','Style','checkbox');
pard.blinkcoded.position=[5,1];
pard.blinkcoded.Width=2;

pard.plugininfo.type='ROI_Analyze';

end