function module=makeSMAPplugin(obj,pluginpath)
%makes a plugin window and initializes it.
module=plugin(pluginpath{:});
name=pluginpath{end};
p.Vrim=100;

module.handle=figure('MenuBar','none','Toolbar','none','Name',name);
module.attachPar(obj.P);
module.attachLocData(obj.locData);
p.Xrim=10;
module.setGuiAppearence(p)
module.makeGui;
end