% calcMAE() - Calculate Mean Average Error in a dataset pre/post- any
%             particular processing step.
%
% Usage: 
%   >> MAE = calcMAE(preEEG, postEEG, order)
%
% Inputs:
%   preEEG  - EEG signal pre-processing step in channels x samples format
%   postEEG - EEG signal post-processing step in channels x samples format
%   order   - The order in which to subtract the pre and post EEGs.
%             {For 1, preEEG - postEEG; For 2, postEEG - preEEG} 
%
% Outputs:
%   MAE - Mean Average Error across all channels and timepoints
%
% Author: Alexa D. Monachino, PINE Lab at Northeastern University, 2021
%         Laurel J. Gabard-Durnam, PINE Lab at Northeastern University, 2021
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

function MAE = calcMAE(preEEG, postEEG, order)
    MAE = zeros(size(preEEG, 1), 1) ;
    for j = 1:size(preEEG, 1)
        mae = 0 ;
        for i = 1:length(preEEG)
            if order == 2; mae = mae + abs(postEEG(j,i) - preEEG(j,i)) ;
            elseif order == 1; mae = mae + abs(preEEG(j,i) - postEEG(j,i)) ;
            end
        end
        [mae] = mae/length(preEEG) ;
        MAE(j) = [mae] ;
    end
    MAE = mean(MAE) ;
end