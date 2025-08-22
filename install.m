function install(savePermanently)
% INSTALL Adds all folders and subfolders in the same directory
% as this file to the MATLAB path.
%
% Usage:
%   install()          % Adds paths for the current MATLAB session only.
%   install(true)      % Adds paths and saves them permanently.
%
% Place this file in your project's root directory.
%
% Input:
%   savePermanently (optional, boolean): If true, saves path for future sessions.
%                                        Defaults to false.

    % Set default for savePermanently if not provided
    if nargin < 1
        savePermanently = false;
    end

    try
        % Get the full path of the current script
        fullPath = mfilename('fullpath');
        
        % Get the parent directory (project's base directory)
        baseDir = fileparts(fullPath);
        
        % Generate path string including base and all subfolders
        pathToAdd = genpath(baseDir);
        
        % Add to MATLAB path
        addpath(pathToAdd);
        
        fprintf('Added "%s" and its subfolders to MATLAB path.\n', baseDir);
        
        % Optionally save path permanently
        if savePermanently
            savepath;
            fprintf('MATLAB path saved permanently.\n');
        else
            fprintf('Path changes are temporary for this session.\n');
            fprintf('To save permanently, call "install(true)" or "savepath" manually.\n');
        end
        
    catch ME
        % Handle errors
        fprintf(2, 'Error during path installation: %s\n', ME.message);
        disp('Check write permissions or file placement.');
    end
end
