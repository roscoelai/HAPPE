% determ_saveFormat() - A helper function for setParams.m from HAPPE. Sets
%                       the save format for the processed files based on
%                       user input through the command line. Options
%                       include .mat, .set, .and .txt. Does not throw
%                       and error on or accept invalid inputs.
%
% Usage: 
%   >> save_as_format = determ_saveFormat(datafileformat, low_density)
%
% Inputs:
%   datafileformat - An integer (0, 1, or 2) representing the original file
%                    format of the data. Output from determ_fileFormat.m
%   low_density    - An integer (0 or 1) representing whether the data is
%                    low-density (1) or high-density (0).
%
% Outputs:
%   save_as_format - An integer (1, 2, or 3) representing the file format
%                    in which to save the processed data.
%
% Author: Alexa D. Monachino, PINE Lab at Northeastern University, 2021
%
% This file is part of HAPPE.
%
% HAPPE is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% HAPPE is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with HAPPE. If not, see <https://www.gnu.org/licenses/>.

function save_as_format = determ_saveFormat(datafileformat, low_density)
    % If the data is in .mat format and low density, the data must be
    % exported in .mat format because channel locations are not included in
    % these files and would not produce proper .set files.
    if datafileformat == 1 && low_density
        disp("Your data will be saved as a .mat files (matlab format)") ;
        save_as_format = 2 ;
    % For all other scenarios, the user can choose .txt, .mat, or .set,
    % with .txt being recommended for ERP timeseries.
    else
        disp("Format to save processed data:") ;
        disp("  1 = .txt file (electrodes as columns, time as rows) - Choose this for ERP timeseries")
        disp("  2 = .mat file (matlab format)") ;
        disp("  3 = .set file (EEGLab format)") ;
        while true
            save_as_format = input('> ') ;
            if ismember(save_as_format, [1, 2, 3]); break ;
            else; disp("Invalid input: please enter 1, 2, or 3.") ;
            end
        end
    end
end