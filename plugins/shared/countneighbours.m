function neighboursout=countneighbours(fp,dx)
if ~isfield(fp,'x')
    x=fp(:,2);y=fp(:,3);
else
    x=fp.x; y=fp.y;
end

[~,ind]=sort(x);

x=x(ind);y=y(ind);

% sdfg
%sort back:
nums=1:length(x);

[~,ind2]=sort(nums(ind));


 neighbours=countneighbours2Dcirc(double(x),double(y),double(dx));
%neighbours=evalc('countneighbours2Dcirc(double(x),double(y),double(dx))');

neighboursout=neighbours(ind2);
%c: n=f(x,y,dx). x already sorted

% indx=1;
% for this=1:10000%length(x)
%     while (x(indx)<(x(this)-dx))&&indx<length(x) %search first x in range
%         indx=indx+1;
%     end
%     testind=indx;
%     while (x(testind)<x(this)+dx)&&testind<length(x) %test all in allowed range
%         if (y(testind)-y(this))^2<dx^2
%             neighbours(this)=neighbours(this)+1;
%         end
%         testind=testind+1;
%     end
% end
