function list=mystruct2list(st,prefix)
if nargin<2
    prefix='';
end
list={};
fn=fieldnames(st);
for k=1:length(fn)
    sh=st.(fn{k});
    if isstruct(sh)
        listnew=mystruct2list(sh,[prefix  fn{k} '-']);
        list=vertcat(list,listnew);
    else
        list(end+1,:)={[prefix fn{k}],sh};
    end
end
        

end