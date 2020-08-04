classdef sendDataToVR<interfaces.DialogProcessor
    % 
    methods
        function obj=sendDataToVR(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_layerson','linewidth_roi','znm_min','znm_max','sr_pixrec','layernames'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)     
            locD = obj.locData;
            layerson=find(obj.getPar('sr_layerson'));
            [locs,~, hroi]=locD.getloc({'xnm','ynm','znm','locprecnm','locprecznm','frame','xnmline','ynmline'},'layer',1,'position','roi');
            % You can also do filtering based on the layers
            %tcpipClient = tcpclient('127.0.0.1',5555,'Timeout',30);
            %request = uint8('ID');
            %write(tcpipClient,request);
            %clear tcpipClient;
            
            %tcpipClient = tcpclient('127.0.0.1',5555,'Timeout',30);
            %data = uint8('ROI');
            %write(tcpipClient,data);
            %clear tcpipClient;
              
            color_field = 0;
            
            switch p.colorfield.selection
                case 'X'
                    color_field = 0; 
                case 'Y'
                    color_field = 1;
                case 'Z'
                    color_field = 2;
                case 'Precision'
                    color_field = 3;
                case 'Z Precision'
                    color_field = 4;
                case 'Frames'
                    color_field = 5;
            end
            
            disp(color_field)
            
            nb_bytes = 65336/(4*numel(fieldnames(locs)))
            
            color_map = locD.P.par.layer1_.content.lut.selection;
            
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
            
            
           
                  
            
            
            
            out.locsData = locs;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef
    % You can implement you gui options here.
    pard.tlt.object=struct('String','Select color coding field','Style','text');
    pard.tlt.position=[1,4];
    
    pard.colorfield.object=struct('String',{{'X','Y','Z','Precision', 'Z Precision','Frames'}},'Style','popupmenu');
    pard.colorfield.position=[2,4];
end