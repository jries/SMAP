classdef (ConstructOnLoad) mParsArgModifiedData < event.EventData
   properties
      modelID
   end
   
   methods
      function data = mParsArgModifiedData(ID)
         data.modelID = ID;
      end
   end
end