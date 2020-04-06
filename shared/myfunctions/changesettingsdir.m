function out = changesettingsdir(file,settingsdir)
if strcmp(file(1:8),'settings')
    out=[settingsdir file(9:end)];
end