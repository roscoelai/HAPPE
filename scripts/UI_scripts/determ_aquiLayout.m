% determ_aquiLayout() - A helper function for setParams.m from HAPPE. Sets
%                       the aquisition layout for the run using user input
%                       from the command line. Does not error out or accept
%                       invalid inputs.
%
% Usage: 
%   >> layout_type = determ_aquiLayout()
%
% Inputs:
%
% Outputs:
%   layout_type - An array representing the layout using the net type, 
%                 represented as a number 1-5, and the number of channels.
%                 [layout-net, number-of-channels]
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

function layout_type = determ_aquiLayout()
    layout_type = zeros(1,2) ;
    while true
        % Prompt user for layout and store input
        layout_type(1,1) = input(['Aquisition layout net type:\n' ...
            '  1 = EGI (Hydrocel) Geodesic Sensor Net\n' ...
            '  2 = Biosemi\n' ...
            '  3 = Brain Products Standard BrainCap (BC)\n' ...
            '  4 = Brain Products Wet-Sponge R-Net for actiCHamp Plus (RNP-AC)\n' ...
            '  5 = Neuroscan Quik-Cap\n' ...
            'NOTE: For other nets, run through BEAPP, as described in the HAPPE manual\n' ...
            '> ']) ;
        % If the entered layout is a valid net, display and ask for the
        % number of channels in the net
        if ismember(layout_type(1,1), [1:5])
            while true
                disp('Number of leads/channels:') ;
                if layout_type(1,1) == 1; fprintf(['For EGI nets, HAPPE supports ' ...
                        '32, 64, 128, and 256 channels.\n']) ;
                elseif layout_type(1,1) == 2; fprintf(['For BioSemi nets, ' ...
                        'HAPPE supports 16, 32, 64, and 128 channels.\n']) ;
                elseif layout_type(1,1) == 3; fprintf(['For BrainProducts ' ...
                        'BC nets, HAPPE supports 22, 32, 64, 96, and 128 ' ...
                        'channels.\n']) ;
                elseif layout_type(1,1) == 4; fprintf(['For BrainProducts ' ...
                        'RNP-AC nets, HAPPE supports 32, 64, 96, and 128 ' ...
                        'channels.\n']) ;
                elseif layout_type(1,1) == 5; fprintf(['For Neuroscan Quik-Caps, ' ...
                        'HAPPE supports 32, 64, and 128 channels.\n']) ;
                end
                layout_type(1,2) = input('> ') ;
                % Confirm that the number of channels selected is valid
                % with the selected layout.
                if (layout_type(1,1) == 1 && ismember(layout_type(1,2), ...
                        [32, 64, 128, 256])) || (layout_type(1,1) == 2 && ...
                        ismember(layout_type(1,2), [16, 32, 64, 128])) || ...
                        (layout_type(1,1) == 3 && ismember(layout_type(1,2), ...
                        [22, 32, 64, 96, 128])) || (layout_type(1,1) == 4 && ...
                        ismember(layout_type(1,2), [32, 64, 96, 128])) || ...
                        (layout_type(1,1) == 5 && ismember(layout_type(1,2), ...
                        [32, 64, 128, 256]))
                    break ;
                else; disp('The entered number of leads is not supported with this net type.') ;
                end
            end
            break ;
        % Otherwise, alert user to invalid input
        else; disp("Invalid input: please enter 1, 2, 3, 4, or 5. Otherwise, please use BEAPP.") ;
        end
    end
end