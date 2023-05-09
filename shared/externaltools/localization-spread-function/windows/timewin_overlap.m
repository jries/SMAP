function overlap = timewin_overlap(tw1, tw2)
i = 1;
j = 1;
sum = 0;

% jump to sane start point
if tw1(i,2) <= tw2(j,1) % this section of tw1 entirely precedes this part of tw2
    i = find(tw1(i:end,2) > tw2(j,1), 1); % find first section of 1 where there is overlap
elseif tw1(i,1) >= tw2(j,2)
    j = find(tw1(1,2) > tw2(:,1), 1); % find first section of 2 where there is overlap
end

while i <= size(tw1, 1) && j <= size(tw2,1)
    sum = sum + max(0, min(tw1(i,2),tw2(j,2)) - max(tw1(i,1),tw2(j,1)));
    if tw1(i,2) <= tw2(j,2) % done with tw1(i,:)
        i = i+1;
    elseif tw2(j,2) <= tw1(i,2) % done with tw2(j,:)
        j = j+1;
    end
end
overlap = sum;