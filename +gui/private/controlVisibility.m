function controlVisibility(obj,fields,state)
if ~iscell(fields)
    fields={fields};
end
for k=1:length(fields)
    obj.(fields{k}).Visible=state;
end