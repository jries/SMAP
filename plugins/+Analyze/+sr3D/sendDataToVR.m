classdef sendDataToVR<interfaces.DialogProcessor
    % 
    methods
        function obj=sendDataToVR(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_layerson','linewidth_roi','znm_min','znm_max','sr_pixrec','layernames'};
            obj.showresults=false;
        end
        
        function out=run(obj,p)     
            locD = obj.locData;
            layerson=find(obj.getPar('sr_layerson'));
            
            tcpipClient = tcpclient('127.0.0.1',5555,'Timeout',30);
            request = uint8('CH');
            write(tcpipClient,request);
            clear tcpipClient;
            
            tcpipClient = tcpclient('127.0.0.1',5555,'Timeout',30);
            data = uint32(length(layerson));
            write(tcpipClient,data);
            clear tcpipClient;
            
            
            for layer = 1:min(2,length(layerson)) %only maximally two channels
                
                colorfieldselected=p.colorfield.selection;
                layerfield=obj.getPar(['layer' num2str(layerson(layer)) '_renderfield']).selection;
                
                [locs,~, hroi]=locD.getloc({'xnm','ynm','znm','locprecnm','locprecznm','frame',colorfieldselected, layerfield,'xnmline','ynmline'},'layer',layerson(layer),'position','roi');
                % You can also do filtering based on the layers
                %tcpipClient = tcpclient('127.0.0.1',5555,'Timeout',30);
                %request = uint8('ID');
                %write(tcpipClient,request);
                %clear tcpipClient;

                %tcpipClient = tcpclient('127.0.0.1',5555,'Timeout',30);
                %data = uint8('ROI');
                %write(tcpipClient,data);
                %clear tcpipClient;
                
%                 color_field = 0;
% 
       color_field=5;
                if p.tlt
              
%                     color_field = p.colorfield.Value - 2 ;
%                     disp(color_field)
%                     if color_field > 5
%                         color_field = 3;
%                     end
                     cfield=colorfieldselected;
                else 
%                     color_field=5;
                    cfield=layerfield;
                end
%                 
%                 switch cfield %in case existing parameters are selected
%                     case 'xnm'
%                         color_field = 0; 
%                     case 'ynm'
%                         color_field = 1;
%                     case 'znm'
%                         color_field = 2;
%                     case 'locprecnm'
%                         color_field = 3;
%                     case 'frame'
%                         color_field = 4;
%                     case 'locprecznm'
%                         color_field = 5;
%                 end
%                 disp(color_field)

                nb_bytes = 65336/(4*numel(fieldnames(locs)))
                
                if layer == 2
                    color_map = locD.P.par.layer2_.content.lut.selection;
                
                else
                    color_map = locD.P.par.layer1_.content.lut.selection;             
                end
                
                disp(color_map);
                
                tot_count=0;

                nb_packages = ceil(length(locs.xnm)/nb_bytes);

                tcpipClient = tcpclient('127.0.0.1',5555,'Timeout',30);
                request = uint8('PK');
                write(tcpipClient,request);
                clear tcpipClient;

                tcpipClient = tcpclient('127.0.0.1',5555,'Timeout',30);
                data = uint32(nb_packages);
                write(tcpipClient,data);
                clear tcpipClient;
                disp(length(locs.xnm))

                for j = 1:nb_packages

                    tcpipClient = tcpclient('127.0.0.1',5555,'Timeout',30);
                    request = uint8('LC');
                    write(tcpipClient,request);
                    clear tcpipClient;

                    package = [];
                    for i=1:nb_bytes

                        tot_count = tot_count + 1;
                        package = [package; locs.xnm(tot_count)];
                        package = [package; locs.ynm(tot_count)];
                        package = [package; locs.znm(tot_count)];
                        package = [package; locs.locprecnm(tot_count)];
                        package = [package; locs.frame(tot_count)];
%                         package = [package; locs.locprecznm(tot_count)]; %this was set to zero?
                        package = [package; single(locs.(cfield)(tot_count))];;%I added two fields
                        %package = [package; locs.(layerfield)(tot_count)]; %this as well
                        
                        %PLACEHOLDER COLUMN
                        package = [package; 0];
                        package = [package; 0];
                        if tot_count == length(locs.xnm)
                            break
                        end 
                    end
                    tcpipClient = tcpclient('127.0.0.1',5555,'Timeout',30);
                    data = single(package);
                    write(tcpipClient, data);
                    clear tcpipClient;

                    if tot_count == length(locs.xnm)
                        break
                    end
                end 

                tcpipClient = tcpclient('127.0.0.1',5555,'Timeout',30);
                request = uint8('CF');
                write(tcpipClient,request);
                clear tcpipClient;

                tcpipClient = tcpclient('127.0.0.1',5555,'Timeout',30);
                data = uint8(color_field);
                write(tcpipClient,data);
                clear tcpipClient;




                tcpipClient = tcpclient('127.0.0.1',5555,'Timeout',30);
                request = uint8('CO');
                write(tcpipClient,request);
                clear tcpipClient;

                tcpipClient = tcpclient('127.0.0.1',5555,'Timeout',30);
                data = uint8(color_map);
                write(tcpipClient,data);
                clear tcpipClient;
             
            end



            out.locsData = locs;
        end
        
            function pard=guidef(obj)
                pard=guidef;
            end
   
    end
end



function pard=guidef
    % You can implement you gui options here.
    
    pard.text.object=struct('String','Send localizations to Genuage. Consult Info for instructions.','Style','text');
    pard.text.position=[1,1];
    pard.text.Width=4;
%         pard.intensitycoding.object=struct('String','Use intensity color coding:','Style','checkbox','Value',1);
%     pard.intensitycoding.position=[2,1];
%     pard.intensitycoding.Width=2;    
    pard.tlt.object=struct('String','Overwrite color coding with:','Style','checkbox');
    pard.tlt.position=[3,1];
    pard.tlt.Width=2;    
    pard.colorfield.object=struct('String',{{'X','Y','Z','Precision', 'Z Precision','Frames'}},'Style','popupmenu');
    pard.colorfield.position=[3,3];
    
    pard.syncParameters={{'locFields','colorfield',{'String'},{}}};
    pard.plugininfo.type='ProcessorPlugin';
    pard.plugininfo.description='Genuage';
end
    

