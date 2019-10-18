classdef CorrectUnconnectbyConnect<interfaces.DialogProcessor
%     Writes the coordiates of the connected (merged) localizations to the
%     correspondng non-connected localizations.
    methods
        function obj=CorrectUnconnectbyConnect(varargin)     
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p) 
            out=[];
            lu=obj.locData.loc;
            lg=obj.locData.grouploc;
            
            [gisort,inds]=sort(lu.groupindex);
%             idx1=1;
            idx2=1;
            for k=1:length(lg.groupindex);
                while idx2<length(gisort)&&gisort(idx2)<lg.groupindex(k)
                    idx2=idx2+1;
                end
                while idx2<length(gisort)&&gisort(idx2)==lg.groupindex(k)
                    obj.locData.loc.xnm(inds(idx2))=lg.xnm(k);
                    obj.locData.loc.ynm(inds(idx2))=lg.ynm(k);
                    if isfield(obj.locData.loc,'znm')
                        obj.locData.loc.znm(inds(idx2))=lg.znm(k);
                    end
                    if isfield(obj.locData.loc,'locprecznm')
                        obj.locData.loc.locprecznm(inds(idx2))=lg.locprecznm(k);
                    end
                    obj.locData.loc.locprecnm(inds(idx2))=lg.locprecnm(k);
                    idx2=idx2+1;
                end
            end
                
             
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Writes the coordiates of the connected (merged) localizations to the correspondng non-connected localizations.';
end