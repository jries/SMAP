% PROCESSMANAGER - Launch and manage external processes
%
%     obj = processManager(varargin);
%
%     Class for launching and managing processes than run asynchronously
%     and in parallel to the main Matlab process. This could be done with 
%     something like 
%     
%     >> system('dir &');
%
%     but using processManager allows you to start and stop processes, peek
%     and check on the progress of running processes, all while allowing you 
%     to continue working in the main Matlab process.
%
%     All inputs are passed in using name/value pairs. The name is a string
%     followed by the value (described below).
%     The only required input is the command.
%     The order of the pairs does not matter, nor does the case.
%
%     More information and can be found on GitHub:
%     https://github.com/brian-lau/MatlabProcessManager/wiki
%
% INPUTS
%     command      - command to execute in separate process, can take the
%                    form of
%                    1) string defining complete command including arguments
%                    2) cell array of strings, parsing the command and each
%                    argument into a separate cell array element
%
% OPTIONAL
%     id           - string identifier for process, default ''
%     workingDir   - string defining working directory
%     envp         - not working yet
%     printStdout  - boolean to print stdout stream, default true
%     printStderr  - boolean to print stderr stream, default true
%     wrap         - number of columns for wrapping lines, default = 80
%     keepStdout   - boolean to keep stdout stream, default false
%     keepStderr   - boolean to keep stderr stream, default false
%     verbose      - boolean to print processManager info, default false
%     autoStart    - boolean to start process immediately, default true
%     pollInterval - double defining polling interval in sec, default 0.5
%                    Take care with this variable, if set too long, one risks
%                    permanently blocking Matlab when streams buffers are not 
%                    drained fast enough. If you don't want to see output, 
%                    it's safer to set printStdout and printStderr false
%
% METHODS
%     start        - start process(es)
%     stop         - stop process(es)
%     check        - check running process(es)
%     block        - block until done
%
% EXAMPLES
%     % 1) Running a simple command
%     p = processManager('command','nslookup www.google.com');
%
%     % 2) Command with ongoing output
%     p = processManager('command','ping www.google.com');
%     % To keep the process running silently,
%     p.printStdout = false;
%     % ... Check back later
%     p.printStdout = true;
%     % Terminate
%     p.stop();
%
%     % 3) Multiples processes
%     p(1) = processManager('id','google','command','ping www.google.com','autoStart',false);
%     p(2) = processManager('id','yahoo','command','ping www.yahoo.com','autoStart',false);
%     p.start()
%     % Tired of hearing about second process
%     p(2).printStdout = false;
%     % ... if you want to hear back later,
%     p(2).printStdout = true;
%     p.stop();
% 
%     $ Copyright (C) 2017 Brian Lau http://www.subcortex.net/ $
%     Released under the BSD license. The license and most recent version
%     of the code can be found on GitHub:
%     https://github.com/brian-lau/MatlabProcessManager

% TODO
% o If streams are stored maybe we need to buffer a finite number of lines
% o Generate unique names for each timer. check using timerfindall, storename
% o cprintf for colored output for each process?

classdef processManager < handle
   properties(SetAccess = public)
      id
      command
      envp
      workingDir

      printStderr
      printStdout
      wrap
      keepStderr
      keepStdout
      autoStart

      verbose
      pollInterval
   end
   properties(SetAccess = private)
      stderr = {};
      stdout = {};
   end
   properties(SetAccess = private, Dependent = true)
      running
      exitValue
   end
   properties(SetAccess = private, Hidden = true)
      process
      state % processState() object
      stderrReader
      stdoutReader
      pollTimer
   end
   properties(SetAccess = protected)
      version = '0.5.2';
   end
   
   methods
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% Constructor
      function self = processManager(varargin)
         % Constructor, arguments are taken as name/value pairs
         %
         % id           - string identifier for process, default ''
         % command      - string defining command to execute, required
         % workingDir   - string defining working directory
         % envp         - not working yet
         % printStdout  - boolean to print stdout stream, default true
         % printStderr  - boolean to print stderr stream, default true
         % wrap         - number of columns for wrapping lines, default = 80
         % keepStdout   - boolean to keep stdout stream, default false
         % keepStderr   - boolean to keep stderr stream, default false
         % autoStart    - boolean to start process immediately, default true
         % verbose      - boolean to print processManager info, default false
         % pollInterval - double defining polling interval in sec, default 0.5
         %                Take care with this variable, if set too long,
         %                runs the risk of blocking Matlab when streams buffers
         %                not drained fast enough
         %                If you don't want to see output, better to set 
         %                printStdout and printStderr false
         %
         p = inputParser;
         p.KeepUnmatched = false;
         p.FunctionName = 'processManager constructor';
         p.addParamValue('id','');
         p.addParamValue('command','');
         p.addParamValue('workingDir','');
         p.addParamValue('envp','');
         p.addParamValue('printStdout',true);
         p.addParamValue('printStderr',true);
         p.addParamValue('keepStdout',false);
         p.addParamValue('keepStderr',false);
         p.addParamValue('wrap',80);
         p.addParamValue('autoStart',true);
         p.addParamValue('verbose',false);
         p.addParamValue('pollInterval',0.05);
         p.parse(varargin{:});
         
         if ~usejava('jvm')
            error([mfilename ' requires Java to run.']);
         end
         
         self.id = p.Results.id;
         self.workingDir = p.Results.workingDir;
         self.envp = p.Results.envp;
         self.printStdout = p.Results.printStdout;
         self.printStderr = p.Results.printStderr;
         self.keepStdout = p.Results.keepStdout;
         self.keepStderr = p.Results.keepStderr;
         self.wrap = p.Results.wrap;
         self.autoStart = p.Results.autoStart;
         self.verbose = p.Results.verbose;
         self.pollInterval = p.Results.pollInterval;

         self.state = processState();
         
         self.command = p.Results.command;
      end
      
      function set.id(self,id)
         if ischar(id)
            self.id = id;
         elseif isscalar(id)
            self.id = num2str(id);
         else
            error('processManager:id:InputFormat','id must be scalar.');
         end
      end
      
      function set.command(self,command)
         if iscell(command)
            % StringTokenizer is used to parse the command based on spaces
            % this may not be what we want, there is an overload of exec()
            % that allows passing in a String array.
            % http://www.mathworks.com/matlabcentral/newsreader/view_thread/308816
            n = length(command);
            cmdArray = javaArray('java.lang.String',n);
            for i = 1:n
               cmdArray(i) = java.lang.String(command{i});
            end
            self.command = cmdArray;
         elseif ischar(command)
            self.command = command;
         elseif isa(command,'java.lang.String[]') || isa(command,'java.lang.String')
            self.command = command;
         else
            error('processManager:command:InputFormat',...
               'command must be a string, cell array of strings, or java.lang.String array.');
         end
         
         if self.autoStart && ~isempty(self.command)
            self.start();
         end
      end
      
      function set.workingDir(self,workingDir)
         if ~ischar(workingDir)
            error('processManager:workingDir:InputFormat',...
               'command must be a string specifying a directory.');
         end
         if isempty(workingDir)
            self.workingDir = pwd;
         elseif exist(workingDir,'dir') == 7
            self.workingDir = workingDir;
         else
            error('processManager:workingDir:InputFormat',...
               'Not a valid directory name.');
         end
      end
      
      function set.envp(self,envp)
         if isempty(envp)
            self.envp = [];
         elseif ischar(envp)
            temp = javaArray('java.lang.String',1);
            temp(1) = java.lang.String(envp);
            self.envp = temp;
         elseif iscell(envp)
            n = length(envp);
            cmdArray = javaArray('java.lang.String',n);
            for i = 1:n
               cmdArray(i) = java.lang.String(envp{i});
            end
            self.envp = cmdArray;
         else
            error('processManager:envp:InputFormat',...
               'command must be a string, cell array of strings, or java.lang.String array.');
         end
      end
      
      function set.printStdout(self,bool)
         if isscalar(bool) && islogical(bool)
            self.printStdout = bool;
            self.updatePollData('printStdout',bool);
         else
            error('processManager:printStdout:InputFormat',...
               'Input must be a scalar logical');
         end
      end
      
      function set.printStderr(self,bool)
         if isscalar(bool) && islogical(bool)
            self.printStderr = bool;
            self.updatePollData('printStderr',bool);
         else
            error('processManager:printStderr:InputFormat',...
               'Input must be a scalar logical');
         end
      end
      
      function set.keepStdout(self,bool)
         if isscalar(bool) && islogical(bool)
            self.keepStdout = bool;
            self.updatePollData('keepStdout',bool);
         else
            error('processManager:keepStdout:InputFormat',...
               'Input must be a scalar logical');
         end
      end
      
      function set.keepStderr(self,bool)
         if isscalar(bool) && islogical(bool)
            self.keepStderr = bool;
            self.updatePollData('keepStderr',bool);
         else
            error('processManager:keepStderr:InputFormat',...
               'Input must be a scalar logical');
         end
      end
      
      function stderr = get.stderr(self)
         if ~isempty(self.state)
            stderr = self.state.stderr;
         else
            stderr = {};
         end
      end
      
      function stdout = get.stdout(self)
         if ~isempty(self.state)
            stdout = self.state.stdout;
         else
            stdout = {};
         end
      end
      
      function set.wrap(self,x)
         if isscalar(x) && isnumeric(x) && (x>0)
            self.wrap = max(1,round(x));
            self.updatePollData('wrap',self.wrap);
         else
            error('processManager:wrap:InputFormat',...
               'Input must be a scalar > 0');
         end
      end
      
      function set.autoStart(self,bool)
         if isscalar(bool) && islogical(bool)
            self.autoStart = bool;
         else
            error('processManager:autoStart:InputFormat',...
               'Input must be a scalar logical');
         end
      end
      
      function set.verbose(self,bool)
         if isscalar(bool) && islogical(bool)
            self.verbose = bool;
            self.updatePollData('verbose',bool);
         else
            error('processManager:verbose:InputFormat',...
               'Input must be a scalar logical');
         end
      end
      
      function set.pollInterval(self,x)
         if isscalar(x) && isnumeric(x) && (x>0)
            self.pollInterval = x;
         else
            error('processManager:pollInterval:InputFormat',...
               'Input must be a scalar > 0');
         end
      end
      
      function updatePollData(self,name,arg)
         for i = 1:numel(self)
            if ~isempty(self(i).pollTimer) && isvalid(self(i).pollTimer)
               dat = get(self.pollTimer,'UserData');
               dat.(name) = arg;
               set(self.pollTimer,'UserData',dat);
            end
         end
      end

      function start(self)
         runtime = java.lang.Runtime.getRuntime();
         for i = 1:numel(self)
            if isempty(self(i).command)
               continue;
            end
            try
               self(i).process = runtime.exec(self(i).command,...
                  self(i).envp,...
                  java.io.File(self(i).workingDir));
               
               % Process will block if streams not drained
               self(i).stdoutReader = java.io.BufferedReader(...
                  java.io.InputStreamReader(self(i).process.getInputStream()));
               self(i).stderrReader = java.io.BufferedReader(...
                  java.io.InputStreamReader(self(i).process.getErrorStream()));

               % Install timer to periodically drain streams
               % http://stackoverflow.com/questions/8595748/java-runtime-exec               
               self(i).pollTimer = timer(...
                  'ExecutionMode','FixedRate',...
                  'BusyMode','queue',...
                  'Period',self(i).pollInterval,...
                  'Name',[self(i).id '-processManager-pollTimer'],...
                  'StartFcn',@processManager.pollTimerStart,...
                  'StopFcn',@processManager.pollTimerStop,...
                  'TimerFcn',@processManager.poll);

               % Load data for the timer to avoid self reference
               pollData.verbose = self(i).verbose;
               pollData.printStderr = self(i).printStderr;
               pollData.printStdout = self(i).printStdout;
               pollData.wrap = self(i).wrap;
               pollData.stderr = self(i).stderr;
               pollData.stdout = self(i).stdout;
               pollData.keepStderr = self(i).keepStderr;
               pollData.keepStdout = self(i).keepStdout;
               pollData.stderrReader = self(i).stderrReader;
               pollData.stdoutReader = self(i).stdoutReader;
               pollData.process = self(i).process;
               pollData.state = self(i).state;
               pollData.state.id = self(i).id;
               set(self(i).pollTimer,'UserData',pollData);
               
               start(self(i).pollTimer);
            catch err
               if any(strfind(err.message,'java.io.IOException: error=2, No such file or directory'))
                  error('processManager:start:InputFormat',...
                     'Looks like command doesn''t exist. Check spelling or path?');
               else
                  rethrow(err);
               end
            end
         end
      end

      function stop(self)
         for i = 1:numel(self)
            if ~isempty(self(i).process)
               self(i).process.destroy();
               self(i).process.getErrorStream().close()
               self(i).process.getInputStream().close()
               self(i).process.getOutputStream().close();
               self(i).stdoutReader.close();
               self(i).stderrReader.close();
            end
         end
      end

      function running = get.running(self)
         if isempty(self.process)
            running = false;
         else
            running = self.isRunning(self.process);
         end
      end
      
      function exitValue = get.exitValue(self)
         if isempty(self.process)
            exitValue = NaN;
         else
            [~,exitValue] = self.isRunning(self.process);
         end
      end
      
      function check(self)
         for i = 1:numel(self)
            if ~self(i).running && isa(self(i).process,'java.lang.Process')
               fprintf('Process %s finished with exit value %g.\n',self(i).id,self(i).exitValue);
            elseif self(i).running && isa(self(i).process,'java.lang.Process')
               fprintf('Process %s is still running.\n',self(i).id);
            else
               fprintf('Process %s has not been started yet.\n',self(i).id);
            end
         end
      end
      
      function block(self,t)
         % Avoid p.waitfor since it hangs with enough output, and
         % unfortunately, Matlab's waitfor does not work?
         % http://undocumentedmatlab.com/blog/waiting-for-asynchronous-events/
         if nargin < 2
            t = self.pollInterval;
         end
         while any([self.running])
            % Matlab pause() has a memory leak
            % http://undocumentedmatlab.com/blog/pause-for-the-better/
            % http://matlabideas.wordpress.com/2013/05/18/take-a-break-have-a-pause/
            java.lang.Thread.sleep(t*1000);
         end
         self.stop();
      end
      
      function delete(self)
         if ~isempty(self.process)
            self.process.destroy();
            self.process.getErrorStream().close()
            self.process.getInputStream().close()
            self.process.getOutputStream().close();
            self.stdoutReader.close();
            self.stderrReader.close();
         end
         if ~isempty(self.pollTimer)
            if isvalid(self.pollTimer)
               stop(self.pollTimer);
            end
         end
         if self.verbose
            fprintf('processManager cleaned up.\n');
         end
      end
   end
   
   methods(Static)
      function pollTimerStart(t,e)
         pollData = get(t,'UserData');
         if pollData.verbose
            fprintf('processManager starting timer for process %s.\n',pollData.state.id)
         end
      end
      
      function pollTimerStop(t,e)
         pollData = get(t,'UserData');
         if pollData.verbose
            fprintf('processManager stopped, uninstalling timer for process %s.\n',pollData.state.id)
         end
         % Update pollData
         [running,exitValue] = processManager.isRunning(pollData.process);
         pollData.state.running = running;
         pollData.state.exitValue = exitValue;
         set(t,'UserData',[]);
         delete(t);
      end
      
      function poll(t,e)
%           return
         pollData = get(t,'UserData');
         try
            stderr = processManager.readStream(pollData.stderrReader);
%             stdout='';
             stdout = processManager.readStream(pollData.stdoutReader);
         catch err
            if any(strfind(err.message,'java.io.IOException: Stream closed'))
               % pass, this can happen when processManager object is
               % stopped or cleared, and should be harmless
               if pollData.verbose
                  fprintf('projectManager timer is polling a closed stream!\n');
               end
            else
               rethrow(err);
            end
         end

         if pollData.printStderr && exist('stderr','var')
            stderr = processManager.printStream(stderr,pollData.state.id,pollData.wrap);
         end
         if pollData.printStdout && exist('stdout','var')
            stdout = processManager.printStream(stdout,pollData.state.id,pollData.wrap);
         end
         
         if pollData.keepStderr && exist('stderr','var')
            pollData.state.updateStderr(stderr);
         end
         if pollData.keepStdout && exist('stdout','var')
            pollData.state.updateStdout(stdout);
         end
         
         % Run timer StopFcn callback if process is done
         running = processManager.isRunning(pollData.process);
         if ~running
            stop(t);
         end
      end
      
      function lines = readStream(stream)
         % This is potentially fragile since ready() only checks whether
         % there is an element in the buffer, not a complete line.
         % Therefore, readLine() can block if the process doesn't terminate
         % all output with a carriage return...
         %
         % Alternatives inlcude:
         % 1) Implementing own low level read() and readLine()
         % 2) perhaps java.nio non-blocking methods
         % 3) Custom java class for spawning threads to manage streams
         lines = {};
         while true
            if stream.ready()
               line = stream.readLine();
               if isnumeric(line) && isempty(line)
                  % java null is empty double in matlab
                  % http://www.mathworks.com/help/matlab/matlab_external/passing-data-to-a-java-method.html
                  break;
               end
               c = char(line);
%                disp(c)
               lines = cat(1,lines,c);
               break
            else
               break;
            end
         end
      end
      
      function str = printStream(c,prefix,wrap)
         if nargin < 2
            prefix = '';
         end
         if nargin < 3
            wrap = 80;
         end
         str = {};
         for i = 1:length(c)
            if exist('linewrap','file') == 2
               if isempty(prefix)
                  tempStr = linewrap(c{i},wrap);
               else
                  tempStr = linewrap([prefix ': ' c{i}],wrap);
               end
            else
               if isempty(prefix)
                  tempStr = c{i};
               else
                  tempStr = [prefix ': ' c{i}];
               end
            end
            str = cat(1,str,tempStr);
         end
         fprintf('%s\n',str{:});
      end
      
      function [bool,exitValue] = isRunning(process)
         try
            exitValue = process.exitValue();
            bool = false;
         catch err
            if any(...
                  [strfind(err.message,'java.lang.IllegalThreadStateException: process hasn''t exited')...
                   strfind(err.message,'java.lang.IllegalThreadStateException: process has not exited')]...
                  )
               bool = true;
               exitValue = NaN;
            else
               rethrow(err);
            end
         end
      end
   end
end
