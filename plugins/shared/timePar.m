classdef timePar
properties
    mainData
end
methods
    function obj = timePar()
       obj.mainData = table([],[],[],[],[],[],[],'VariableNames',{'time', 'th', 'de', 'zm', 'shi', 'nm', 'sha'});
    end
    function plot(obj)

    end
end
end