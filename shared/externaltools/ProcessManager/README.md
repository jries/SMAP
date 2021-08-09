# processManager

A Matlab class for launching and managing processes that run asynchronously from the main Matlab process. This can already be done with something like `system('dir &');` but processManager makes it easy to:

* launch and manage multiple processes
* check on the progress of running processes
* capture & display `stdout` and `stderr` streams of each process
* issue event notifications when processes finish

while allowing you to continue working in the main Matlab process.

Some toy examples are illustrated below and in the [wiki](https://github.com/brian-lau/MatlabProcessManager/wiki). A more elaborate application is a [Matlab interface](https://github.com/brian-lau/MatlabStan) that does MCMC sampling using [Stan](http://mc-stan.org/).

## Installation & Examples
Download [processManager](https://github.com/brian-lau/MatlabProcessManager/archive/master.zip), add the m-file to your Matlab path, and you're ready to go.

processManager was developed and tested on OSX with Matlab 2012a, but should work on all platforms that Matlab supports, so long as it is running >=R2008a (for handle objects) with JDK >=1.1 (this will always be true unless you changed the default JDK).

### Optional
Installing Steve Eddins's [linewrap](http://www.mathworks.com/matlabcentral/fileexchange/9909-line-wrap-a-string) function is useful for dealing with unwrapped messages. His [xUnit test framework](http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework) is required if you want to run the unit tests.

### Examples

#### Running a simple command
```
p = processManager('command','nslookup www.google.com');
```

#### Command with ongoing output
```
p = processManager('command','ping www.google.com');

% To keep the process running silently,
p.printStdout = false;

% ... Check process status
p.check();

% When you want to see the io stream again
p.printStdout = true;

% Terminate
p.stop();
```

#### Multiples processes using object arrays
You can pack multiple processes into an object array for easy management.
```
p(1) = processManager('id','google','command','ping www.google.com','autoStart',false);
p(2) = processManager('id','yahoo','command','ping www.yahoo.com','autoStart',false);
p.start();

% Tired of hearing about second process
p(2).printStdout = false;

% You can check the status of individual processes
p(2).check();

% Or all processes
p.check()

% ... if you want to hear about an individual process later,
p(2).printStdout = true;

% Terminate all processes
p.stop();
```

## Need help?
You may be able to find a solution in the [wiki](https://github.com/brian-lau/MatlabProcessManager/wiki/Potential-gotchas). Otherwise, open an [issue](https://github.com/brian-lau/MatlabProcessManager/issues).

Contributions
--------------------------------
Copyright (c) 2017 Brian Lau [brian.lau@upmc.fr](mailto:brian.lau@upmc.fr), see [LICENSE](https://github.com/brian-lau/MatlabProcessManager/blob/master/LICENSE.txt)

Please feel free to [fork](https://github.com/brian-lau/MatlabProcessManager/fork) and contribute!
