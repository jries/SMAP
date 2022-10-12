function out=mystructdialog(varargin)
in=varargin{1};
if nargin<2
    title='Parameters';
else
    title=varargin{2};
end
out=[];
fn=fieldnames(in);
for k=length(fn):-1:1
    fields{k}=fn{k};
    defAns{k}=converttostring(in.(fn{k}));
end
answer=inputdlg(fields,title,1,defAns);
if ~isempty(answer)
    for k=1:length(fn)
        if isnumeric(in.(fn{k}))||islogical(in.(fn{k}))
            out.(fn{k})=str2num(answer{k});
        else
            out.(fn{k})=(answer{k});
        end
    end
end

end

function out=converttostring(in)
if iscell(in)
    out=join(in,',');
    out=out{1};
elseif ischar(in)
    out=in;
else
    out=num2str(in(:)');
end
end