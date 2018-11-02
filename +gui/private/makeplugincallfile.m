function d=makeplugincallfile(plugindir)
if ~isdeployed
    addpath('plugins')
    addpath('plugins/shared');
%     addpath('settings/workflows');
    
    try
    pold=plugin;
    catch
        pold=[];
    end
%     outputdir='+plugintemp';
%     if ~exist(outputdir,'dir')
%         mkdir(outputdir)
%         rehash
%     end
    
    strcase={};
    strlist={};
    d=dirrecursive(plugindir,2);
    fn1=fieldnames(d);
    for k1=1:length(fn1)
        fn1h=fn1{k1};
        fn2=fieldnames(d.(fn1h));
        for k2=1:length(fn2)
            fn2h=fn2{k2};
            dirh=[plugindir filesep '+' fn1h filesep '+' fn2h filesep];
            fn3=dir([dirh '*.m']);
%             dirh=d.(fn1h).(fn2h);
            if ~isempty(fn3)
                for k3=1:length(fn3)
                    [~,fn3h]=fileparts(fn3(k3).name);
                    ismodule=false;
                    try
                        ptry=pold.(fn1h).(fn2h).(fn3h);
                        ismodule=true;
%                         if length(ptry)>3
                            pname=ptry{4};
%                         else
%                             pname=ptry{3};
%                         end
                        
%                         if length(ptry)>4
                            pkind=ptry{5};
%                         else
%                             pkind='';
%                         end                        
                    catch
                        try
                            module=callmodule(fn1h,fn2h,fn3h);
                            
                        if isa(module,'interfaces.GuiModuleInterface')
                            module.pluginpath={fn1h,fn2h,fn3h};
                            ismodule=true;
                            try
                                guidef=module.guidef;
                                pname=guidef.plugininfo.name;
                            catch
                            pname=module.info.name;
                            end
                            try
                                guidef=module.guidef;
                                pkind=guidef.plugininfo.type;
                            catch
                                pkind='';
                            end
                        end
                        catch err
                            if ~strcmp(err.identifier,'MATLAB:scriptNotAFunction')
                                warning(err.message)
                                disp([fn1h,fn2h,fn3h])
                            end
                        end
                    end
                    if ismodule
                        strcase(end+1:end+2)=makecase(fn1h,fn2h,fn3h);
                        strlist(end+1)=makelist(fn1h,fn2h,fn3h,pname,pkind);
                    end
%                     catch
%                          fn3h
%                     end
                    
                end
            end
%             end
        end
 
    end
    %assemble outher lines
    strpre1{1}='function out=plugin(fn1,fn2,fn3,varargin)';
    strpre1{2}='if nargin>0';
    strpre1{3}='fstr=[fn1 ''.'' fn2 ''.'' fn3];';
    strpre1{4}='switch fstr';
    
    strpost1{1}='end';
    strpost1{2}='module.pluginpath={fn1,fn2,fn3};';
    strpost1{3}='out=module;';
    strpost1{4}='else';
    
    strpost2{1}='end';
    
    fileout='plugins/plugin.m';
    fid=fopen(fileout,'w');
    writecell(fid,strpre1)
    writecell(fid,strcase)
    writecell(fid,strpost1)
    writecell(fid,strlist)
    writecell(fid,strpost2)
    fclose(fid);
     rehash
else
    d='';
end
end

function dirout=dirrecursive(plugindir,depth)
alldir=dir([plugindir filesep '+*']);
    dirout=[];
    if depth>0
        for k=1:length(alldir)
            if alldir(k).isdir && ~strcmp(alldir(k).name(1),'.')
                dirout.(alldir(k).name(2:end))=dirrecursive([plugindir filesep alldir(k).name],depth-1);
            end
        end
    end
end


function out=makecase(a,b,c)
out{1}=['case ''' a '.' b '.' c ''''];
out{2}=['   module=' a '.' b '.' c '(varargin{:});'];
end

function out=makelist(a,b,c,d,e)
if nargin<4
    out{1}=['out.' a '.' b '.' c '={''' a ''',''' b ''',''' c '''};'];
elseif nargin<5
    out{1}=['out.' a '.' b '.' c '={''' a ''',''' b ''',''' c ''',''' d '''};'];
else
    out{1}=['out.' a '.' b '.' c '={''' a ''',''' b ''',''' c ''',''' d ''',''' e '''};'];
end
end

function out=callmodule(a,b,c)
out=eval([ a '.' b '.' c ]);
end

function writecell(fid,str)
for k=1:length(str)
    fprintf(fid,'%s \n',str{k});
end
end