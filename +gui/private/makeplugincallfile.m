function d=makeplugincallfile(plugindir)
if ~isdeployed
    addpath('plugins')
    addpath('plugins/shared');
    
    try
    pold=plugin;
    catch
        pold=[];
    end

    
    strcase={};
    strlist={};
    d=dirrecursive(plugindir,3);
    fn1=fieldnames(d);
    for k1=1:length(fn1)
        fn1h=fn1{k1};
        fn2=fieldnames(d.(fn1h));
        for k2=1:length(fn2)
            fn2h=fn2{k2};
            [strcaseh, strlisth]=makepluginlist(pold,{fn1h,fn2h},plugindir);
            strlist(end+1:end+length(strlisth))=strlisth;
            strcase(end+1:end+length(strcaseh))=strcaseh;
            if ~isempty(d.(fn1h).(fn2h))
                fn3=fieldnames(d.(fn1h).(fn2h));
                for k3=1:length(fn3)
                    fn3h=fn3{k3};
                    [strcaseh, strlisth]=makepluginlist(pold,{fn1h,fn2h,fn3h},plugindir);
                    strlist(end+1:end+length(strlisth))=strlisth;
                    strcase(end+1:end+length(strcaseh))=strcaseh;
                end
                %another iteration
            end
        end
    end
    %assemble outher lines
    strpre1{1}='function out=plugin(fn1,fn2,fn3,fn4,varargin)';
    strpre1{2}='if nargin>0';
    strpre1{3}='fstr=[fn1 ''.'' fn2 ''.'' fn3 ''.'' fn4];';
    strpre1{4}='switch fstr';
    
    strpost1{1}='end';
    strpost1{2}='module.pluginpath={fn1,fn2,fn3,fn4};';
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

function [strcase, strlist]=makepluginlist(pold,fn,plugindir)
strcase={};
strlist={};
dirh=[plugindir filesep '+' fn{1} filesep '+' fn{2} filesep];
if length(fn)>2
    dirh=[dirh '+' fn{3} filesep];
else
    fn{3}='x';
end
fn3=dir([dirh '*.m']);
if ~isempty(fn3)
    for k3=1:length(fn3)
        [~,fn3h]=fileparts(fn3(k3).name);
        ismodule=false;
        try
            ptry=pold.(fn{1}).(fn{2}).(fn{3}).(fn3h);
            ismodule=true;
            pname=ptry{5};
            pkind=ptry{6};                      
        catch
            try
                module=callmodule(fn{1},fn{2},fn{3},fn3h);
                
            if isa(module,'interfaces.GuiModuleInterface')
                module.pluginpath={fn{1},fn{2},fn{3},fn3h};
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
                    disp([fn{1},fn{2},fn{3},fn3h])
                end
            end
        end
        if ismodule
            strcase(end+1:end+2)=makecase(fn{1},fn{2},fn{3},fn3h);
            strlist(end+1)=makelist(fn{1},fn{2},fn{3},fn3h,pname,pkind);
        end
%                     catch
%                          fn3h
%                     end
        
    end
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


function out=makecase(a,b,c,d)
out{1}=['case ''' a '.' b '.' c '.' d ''''];
if c=='x'
    out{2}=['   module=' a '.' b '.' d '(varargin{:});'];
else
    out{2}=['   module=' a '.' b '.' c '.' d '(varargin{:});'];
end
end

% function out=makelist(a,b,c,d,e,f)
% if nargin<4
%     out{1}=['out.' a '.' b '.' c '={''' a ''',''' b ''',''' c '''};'];
% elseif nargin<5
%     out{1}=['out.' a '.' b '.' c '={''' a ''',''' b ''',''' c ''',''' d '''};'];
% elseif nargin<6
%     out{1}=['out.' a '.' b '.' c '={''' a ''',''' b ''',''' c ''',''' d ''',''' e '''};'];
% else
%         out{1}=['out.' a '.' b '.' c '={''' a ''',''' b ''',''' c ''',''' d ''',''' e ''',''' f '''};'];
% end
% end
function out=makelist(a,b,c,d,e,f)
out{1}=['out.' a '.' b '.' c '.' d '={''' a ''',''' b ''',''' c ''',''' d ''',''' e ''',''' f '''};'];
end

function out=callmodule(a,b,c,d)
if c=='x'
    out=eval([ a '.' b '.' d ]);
else
    out=eval([ a '.' b '.' c '.' d ]);    
end
end

function writecell(fid,str)
for k=1:length(str)
    fprintf(fid,'%s \n',str{k});
end
end