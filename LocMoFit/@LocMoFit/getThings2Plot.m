function items = getThings2Plot(obj)

for m = obj.numOfModel:-1:1
    mPar = obj.exportPars(m,'mPar');
    oneItems = obj.model{m}.modelObj.getThings2Plot(mPar);
    items{m} = oneItems;
end

end