function uninstall(savePermanently)
% UNINSTALL Removes all folders and subfolders in the same directory
% as this file from the MATLAB path.
%
% Usage:
%   uninstall()          % Removes paths for the current MATLAB session only.
%   uninstall(true)      % Removes paths and saves the change permanently.
%
% Place this file in your project's root directory, next to 'install.m'.
%
% Input:
%   savePermanently (optional, boolean): If true, the path changes will be
%                                        saved for future MATLAB sessions.
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
        
        % Generate the path string that was previously added
        pathToRemove = genpath(baseDir);
        
        % Remove paths from MATLAB path
        rmpath(pathToRemove);
        
        fprintf('Successfully removed "%s" and its subfolders from MATLAB path.\n', baseDir);
        
        % Optionally save path permanently
        if savePermanently
            savepath;
            fprintf('MATLAB path changes saved permanently.\n');
        else
            fprintf('Path changes are temporary for this session.\n');
            fprintf('To save permanently, call "uninstall(true)" or "savepath" manually.\n');
        end
        
    catch ME
        % Handle errors
        fprintf(2, 'Error during path uninstallation: %s\n', ME.message);
        disp('Please ensure you have appropriate write permissions if attempting to save the path permanently.');
        disp('Also, verify that the "uninstall.m" file is correctly placed in your project''s root directory.');
    end
end
