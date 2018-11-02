function choice=nonmodaldialog(title,texta, textb)
choice=[];
h=dialog('Name',title,'WindowStyle','normal');
h.Position(3)=200;h.Position(4)=75;
usebtn=uicontrol('Parent',h,'Position',[2 60 195 20],'String',title,'Style','text');
usebtn=uicontrol('Parent',h,'Position',[2 2 90 30],'String',texta,'Callback',@btcallback);
cbtn=uicontrol('Parent',h,'Position',[102 2 90 30],'String',textb,'Callback',@btcallback);
uiwait(h);
    function btcallback(a,b)
        choice=a.String;
        delete(h);
    end
end