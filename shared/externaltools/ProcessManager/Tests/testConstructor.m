% Requires the xUnit test framework
% http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework

function test_suite = testConstructor
initTestSuite;

function testNoArgs
p = processManager();
assertEqual(p.running,false);
assertEqual(p.exitValue,NaN);
assertTrue(isa(p,'processManager'),'Constructor failed to create pointProcess without inputs');

function testArgs
p = processManager('id',999,...
                   'command','ping www.google.com',...
                   'workingDir',tempdir,...
                   'envp',{'FOO=test'},...
                   'printStdout',false,...
                   'printStderr',false,...
                   'keepStdout',true,...
                   'keepStderr',true,...
                   'wrap',99,...
                   'autoStart',false,...
                   'verbose',true,...
                   'pollInterval',1 ...
                   );
assertEqual(p.id,'999');
assertEqual(p.command,'ping www.google.com');
assertEqual(p.workingDir,tempdir);
%
assertEqual(p.printStdout,false);
assertEqual(p.printStderr,false);
assertEqual(p.keepStdout,true);
assertEqual(p.keepStderr,true);
assertEqual(p.wrap,99);
assertEqual(p.autoStart,false);
assertEqual(p.verbose,true);
assertEqual(p.pollInterval,1);

function testBadArgs
f = @() processManager('id',{'should' 'not' 'work'});
assertExceptionThrown(f,'processManager:id:InputFormat');

f = @() processManager('command',1);
assertExceptionThrown(f, 'processManager:command:InputFormat');
f = @() processManager('command',{1 2});
assertExceptionThrown(f, 'processManager:start:InputFormat');
f = @() processManager('command','thiscommanddoesnotexist');
assertExceptionThrown(f, 'processManager:start:InputFormat');
f = @() processManager('command',{'this' 'command' 'does' 'not' 'exist'});
assertExceptionThrown(f, 'processManager:start:InputFormat');

f = @() processManager('workingDir',1);
assertExceptionThrown(f, 'processManager:workingDir:InputFormat');
f = @() processManager('workingDir','/this/directory/does/not/exist/');
assertExceptionThrown(f, 'processManager:workingDir:InputFormat');

f = @() processManager('printStderr',1);
assertExceptionThrown(f, 'processManager:printStderr:InputFormat');
f = @() processManager('printStdout',1);
assertExceptionThrown(f, 'processManager:printStdout:InputFormat');
f = @() processManager('printStderr','t');
assertExceptionThrown(f, 'processManager:printStderr:InputFormat');
f = @() processManager('printStdout','t');
assertExceptionThrown(f, 'processManager:printStdout:InputFormat');
f = @() processManager('printStderr',{true true});
assertExceptionThrown(f, 'processManager:printStderr:InputFormat');
f = @() processManager('printStdout',{true true});
assertExceptionThrown(f, 'processManager:printStdout:InputFormat');

f = @() processManager('keepStderr',1);
assertExceptionThrown(f, 'processManager:keepStderr:InputFormat');
f = @() processManager('keepStdout',1);
assertExceptionThrown(f, 'processManager:keepStdout:InputFormat');
f = @() processManager('keepStderr','t');
assertExceptionThrown(f, 'processManager:keepStderr:InputFormat');
f = @() processManager('keepStdout','t');
assertExceptionThrown(f, 'processManager:keepStdout:InputFormat');
f = @() processManager('keepStderr',{true true});
assertExceptionThrown(f, 'processManager:keepStderr:InputFormat');
f = @() processManager('keepStdout',{true true});
assertExceptionThrown(f, 'processManager:keepStdout:InputFormat');

f = @() processManager('wrap',0);
assertExceptionThrown(f, 'processManager:wrap:InputFormat');
f = @() processManager('wrap',{100});
assertExceptionThrown(f, 'processManager:wrap:InputFormat');
f = @() processManager('wrap',[10 10]);
assertExceptionThrown(f, 'processManager:wrap:InputFormat');

f = @() processManager('autoStart',0);
assertExceptionThrown(f, 'processManager:autoStart:InputFormat');
f = @() processManager('autoStart','t');
assertExceptionThrown(f, 'processManager:autoStart:InputFormat');
f = @() processManager('autoStart',{0});
assertExceptionThrown(f, 'processManager:autoStart:InputFormat');

f = @() processManager('verbose',0);
assertExceptionThrown(f, 'processManager:verbose:InputFormat');
f = @() processManager('verbose','true');
assertExceptionThrown(f, 'processManager:verbose:InputFormat');
f = @() processManager('verbose',[true true]);
assertExceptionThrown(f, 'processManager:verbose:InputFormat');

f = @() processManager('pollInterval',0);
assertExceptionThrown(f, 'processManager:pollInterval:InputFormat');
f = @() processManager('pollInterval','t');
assertExceptionThrown(f, 'processManager:pollInterval:InputFormat');
f = @() processManager('pollInterval',{0});
assertExceptionThrown(f, 'processManager:pollInterval:InputFormat');

