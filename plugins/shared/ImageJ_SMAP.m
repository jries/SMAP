function ImageJ_SMAP(open_imagej, verbose,imagej_directory_in)
    %% This script adds the ImageJ libraries to the MATLAB classpath. By default,
    %% it will also start up an ImageJ instance within MATLAB.
    %% Parameters:
    %%  open_imagej - If false, an ImageJ instance will not be launched. Default: true
    %%  verbose - If true, a confirmation message will be printed the first time
    %%            a jar is added to the MATLAB classpath. Default: false
    %% Author: Jacques Pecreaux, Johannes Schindelin, Jean-Yves Tinevez, Mark Hiner

    if nargin < 1
        open_imagej = true;
    end

    if nargin < 2
        verbose = false;
    end;

    %% Get the ImageJ directory
%     imagej_directory = fileparts(fileparts(mfilename('fullpath')));

    %% Get the Java classpath
    classpath = javaclasspath('-all');

    %% Add all libraries in jars/ and plugins/ to the classpath

    % Switch off warning
    warning_state = warning('off');
    imagej_directory=fileparts(imagej_directory_in);
    add_to_classpath(classpath, fullfile(imagej_directory,'jars'), verbose);
    add_to_classpath(classpath, fullfile(imagej_directory,'plugins'), verbose);

    % Switch warning back to initial settings
    warning(warning_state)

    % Set the ImageJ directory (and plugins.dir which is not ImageJ.app/plugins/)
    java.lang.System.setProperty('ij.dir', imagej_directory);
    java.lang.System.setProperty('plugins.dir', imagej_directory);

    %% Maybe open the ImageJ window
    import net.imagej.matlab.*;
    if open_imagej
        ImageJMATLAB.start();
    else
        % initialize ImageJ with the headless flag
        ImageJMATLAB.start('--headless');
    end

    % Make sure that the scripts are found.
    % Unfortunately, this causes a nasty bug with MATLAB: calling this
    % static method modifies the static MATLAB java path, which is
    % normally forbidden. The consequences of that are nasty: adding a
    % class to the dynamic class path can be refused, because it would be
    % falsy recorded in the static path. On top of that, the static
    % path is fsck in a weird way, with file separator from Unix, causing a
    % mess on Windows platform.
    % So we give it up as now.
    % %    imagej.User_Plugins.installScripts();
end

function add_to_classpath(classpath, directory, verbose)
    % Get all .jar files in the directory
    dirData = dir(directory);
    dirIndex = [dirData.isdir];
    jarlist = dir(fullfile(directory,'*.jar'));
    path_= cell(0);
    for i = 1:length(jarlist)
      if not_yet_in_classpath(classpath, jarlist(i).name)
        if verbose
          disp(strcat(['Adding: ',jarlist(i).name]));
        end
        path_{length(path_) + 1} = fullfile(directory,jarlist(i).name);
      end
    end

    %% Add them to the classpath
    if ~isempty(path_)
        javaaddpath(path_, '-end');
    end

    %# Recurse over subdirectories
    subDirs = {dirData(dirIndex).name};
    validIndex = ~ismember(subDirs,{'.','..'});

    for iDir = find(validIndex)
      nextDir = fullfile(directory,subDirs{iDir});
      add_to_classpath(classpath, nextDir, verbose);
    end
end

function test = not_yet_in_classpath(classpath, filename)
%% Test whether the library was already imported
expression = strcat([filesep filename '$']);
test = isempty(cell2mat(regexp(classpath, expression)));
end
