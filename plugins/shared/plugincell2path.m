function file = plugincell2path( inputc )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
file=['plugins' filesep '+' inputc{1} filesep '+' inputc{2} filesep inputc{3}];

end

