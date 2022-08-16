function polyout = polylineToMask(inSpots, varargin) %, 
% Makes a polygon around specified polyline with a certain margin

defaultMargin = 0.25;

p = inputParser;
p.addParameter('margin',defaultMargin);
p.parse(varargin{:});

try
    polyout = polybuffer(inSpots,'lines',margin);

catch ME
    warning('Error in transforming the polyline to polygon');
    print(ME);

end
end