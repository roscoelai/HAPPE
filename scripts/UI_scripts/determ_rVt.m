% determ_rVt() - A helper function for setParams.m from HAPPE. Determines
%                whether the paradigm used for the data is
%                resting-state/baseline data or task-related data. In the
%                case of task-related data, will also ask if performing ERP
%                analyses. If resting-state, ERP analysis is automatically
%                disabled.
%
% Usage: 
%   >> [task_EEG_processing, ERP_analysis] = determ_rVt()
%
% Inputs:
%
% Outputs:
%   task_EEG_processing - A binary integer of 1 or 0 that indicates whether
%                         the data is task-related or resting-state,
%                         respectively.
%   ERP_analysis        - A binary interger of 1 or 0 that indicates
%                         whether ERP analyses are being performed or not,
%                         respectively. If so, it activates the HAPPE+ER
%                         pipeline. Is automatically set to 0 in
%                         resting-state data.
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

function [task_EEG_processing, ERP_analysis] = determ_rVt()
    disp("Enter data type:") ;
    disp("  rest = Resting-State EEG") ;
    disp("  task = Task-Related EEG") ;
    task_EEG_processing = choose2("rest", "task") ;
    if task_EEG_processing
        disp("Performing event-related potential (ERP) analysis? [Y/N]") ;
        ERP_analysis = choose2("n", "y") ;
    else; ERP_analysis = 0 ;
    end
end