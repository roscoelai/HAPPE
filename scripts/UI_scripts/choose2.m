% choose2() - A helper function for HAPPE, allowing the user to make a
%             choice between two options via command-line input without
%             throwing an error or incorrectly accepting invalid input.
%             Is not case-sensitive.
%
% Usage: 
%   >> user_answ = choose2(choice1, choice2)
%
% Inputs:
%   choice1 - The first user option in the form of a string/char array.
%             This option, when selected, is coded as a 0.
%   choice2 - The second user option in the form of a string/char array.
%             This option, when selected, is coded as a 1.
%
% Outputs:
%   user_answ - A binary value of 0 or 1, reflecting the user's selection
%               of choice1 or choice2, respectively.
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

function user_answ = choose2(choice1, choice2)
    while true
        user_answ = input('> ', 's') ;
        if strcmpi(user_answ, choice1); user_answ = 0 ; break ;
        elseif strcmpi(user_answ, choice2); user_answ = 1 ; break ;
        else; fprintf("Invalid input: please enter %s or %s\n", choice1, choice2) ;
        end
    end
end