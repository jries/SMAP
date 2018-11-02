function badind=removesingles(positions,radius, minneighbours,iterations)

for k=1:iterations
    nb=countneighbours(positions,radius);
    badind=(nb<=minneighbours);
end