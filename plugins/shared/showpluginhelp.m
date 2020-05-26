function hann=showpluginhelp(hp,td,pos,fs)
if nargin<3||isempty(pos)
    pos=[0 0 1 1];
end
if nargin<4||isempty(fs)
    fs=12;
end

 if isstruct(td) %different interpreters
    fn=fieldnames(td);
    ind=(contains(fn,'latex'));
    interpreter0=fn{~ind};
    if ~any(ind) %only tex or none
%         split=false;
        text0=td.(interpreter0);
    elseif length(fn)==1 %only latex
%         split=false;
        text0=td.latex;
        interpreter0='latex';
    else
%         split=true;

        text0=td.(interpreter0);
        textlatex=td.latex;
        ratio=min(0.75,max(.25,2*length(textlatex)/(length(text0)+length(textlatex))));
        postex=pos;
        postex(4)=ratio;
        pos(2)=ratio;pos(4)=(1-ratio);

        hlatex=annotation(hp,'textbox',postex,...
         'FontSize',fs,'HorizontalAlignment','left',...
         'BackgroundColor','w','FitBoxToText','off','EdgeColor','w',...
         'String',formattext(textlatex),'Interpreter','latex','VerticalAlignment','middle');
        hann{2}=hlatex;
    end
else
    text0=td;
end

  htxt=annotation(hp,'textbox',pos,...
     'FontSize',fs,'HorizontalAlignment','left',...
     'BackgroundColor','w','FitBoxToText','off','EdgeColor','w',...
     'String',formattext(text0),'Interpreter',interpreter0);
  htxt.Position=pos;
 hann{1}=htxt;
end

function txt=formattext(td)
 if ~iscell(td)
  txt=strrep(td,char(9),' ');
  txt=strrep(txt,'\n',newline);
 else
     txt=td;
 end
end