% PROCESSSTATE - Utility class for processManager
% 
% Useful for sending notifications based on exitValue of process.
%
classdef processState < handle
   properties
      id
      running = 0
      exitValue = NaN
      stderr = {}
      stdout = {}
   end
   events
      exit
      exit_success
      exit_failure
   end
   
   methods
      function self = processState(varargin)
      end
      
      function set.running(self,val)
         self.running = val;
      end
      
      function set.exitValue(self,val)
         self.exitValue = val;
         if val >= 0
            notify(self,'exit');
         end
         if val == 0
            notify(self,'exit_success');
         end
         if val > 0
            notify(self,'exit_failure');
         end
      end
      
      function updateStderr(self,msg)
         self.stderr = cat(1,self.stderr,msg);
      end
      
      function updateStdout(self,msg)
         self.stdout = cat(1,self.stdout,msg);
      end
   end
end