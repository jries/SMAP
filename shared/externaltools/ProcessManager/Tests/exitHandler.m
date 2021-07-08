function exitHandler(src,data)
   fprintf('\n');
   fprintf('Listener notified!\n');
   fprintf('Process %s exited with exitValue = %g\n',src.id,src.exitValue);
   fprintf('Event name %s\n',data.EventName);
   fprintf('\n');
end