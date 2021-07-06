% matFormat_Task() - A helper function for setParams.m and HAPPE that
%                    collects the path to where task information is stored
%                    through user input in the command line.                   
%
% Usage: 
%   >> matFormat_Task()
%
% Inputs:
%   params - A struct containing all of the parameters needed to run HAPPE.
%
% Outputs:
%   task_event_loc - The path to the event files containing the event
%                    information.
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

function task_event_loc = matFormat_Task()
    while true
        task_event_loc = input(['Path to .txt files containing task event ' ...
            'info:\n> '], 's') ;
        if exist(task_event_loc, 'dir') == 7; break ;
        else; fprintf(['Invalid input: please enter the correct path to ' ...
                'your task event info.\n']) ; 
        end
    end
end