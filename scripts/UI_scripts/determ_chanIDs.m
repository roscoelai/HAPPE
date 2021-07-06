% determ_chanIDs() - A helper function for setParams.m from HAPPE. Sets
%                    the channels of interest for the run through user
%                    input from the command line. Allows the user to select
%                    all channels, or to choose which subset of channels to
%                    include/exclude. For high-density data, the 10-20
%                    channels are automatically included in the channel
%                    IDs.
%
% Usage: 
%   >> [chans_all, chan_IDs] = determ_chanIDs(low_density)
%
% Inputs:
%   low_density    - An integer (0 or 1) representing whether the data is
%                    low-density (1) or high-density (0).
%
% Outputs:
%   chans_all - String indicating whether the user is using all channels or
%               a subset of channels. If selecting a subset, additionally
%               indicates if the chan_IDs should be included or excluded.
%   chan_IDs  - Cell array, with each cell being the name of a channel of
%               interest. The actual contents will vary based on the input
%               parameters and the chans_all selection to best function
%               with HAPPE's current code.
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

function [chans_all, chan_IDs] = determ_chanIDs(low_density)
disp("Examine all channels (all) or only channels of interest (coi)?") ;
while true
    % collect and store user input
    chans_all = input('> ', 's') ;
    % if the user requests all channels...
    if strcmpi(chans_all, 'all')
        chan_IDs = {} ;
        break ;
   % If the user wants a subset, request user input to collect channels of interest
    elseif strcmpi(chans_all, "coi")
        fprintf(['Choose an option for entering channels:\n  include - ' ...
            'Include ONLY the entered channel names.\n  exclude - Include ' ...
            'every channel EXCEPT the entered channel names.\n']) ;
        while true
            chans_all = ['coi_' input('> ', 's')] ;
            if strcmpi(chans_all, 'coi_include') || strcmpi(chans_all, ...
                    'coi_exclude'); break ;
            else
                fprintf(['Invalid input: please enter "include" or ' ...
                    '"exclude" (without quotations).\n']) ;
            end
        end
        % Collect channel names using user input
        fprintf(['Enter channels, including the preceding letter, one at ' ...
            'a time.\nPress Enter/Return between each entry.\nExamples: ' ...
            'E17\n          M1\nWhen you have entered all channels, input ' ...
            '"done" (without quotations).\n']) ;
        if strcmpi(chans_all, 'coi_include') && ~low_density
            disp("NOTE: 10-20 channels are already included.") ;
            chan_IDs = {'FP1' 'FP2' 'F7' 'F3' 'F4' 'F8' 'C3' 'C4' 'T5' 'PZ' ...
                'T6' 'O1' 'O2' 'T3' 'T4' 'P3' 'P4' 'Fz' 'CZ'} ;
            indx = 20 ;
        else
            chan_IDs = {} ;
            indx = 1 ;
        end
        chan_IDs = unique(userInput_cellArray(indx, chan_IDs), 'stable') ;
        break ;
    else; disp("Invalid input: please enter 'all' or 'coi' (without quotations)") ;    
    end
end