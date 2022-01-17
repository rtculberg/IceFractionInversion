%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Riley Culberg
% Date: 4/16/2021
%
% This function chooses the number and area of a set of firn patches that
% covers the desired fraction of the total area without exceeding it.
% 
% Input Variables:
% f - fractional area to tile (scalar double)
% A_tot - total area (scalar double)
% A_min - minimum area of patches (scalar double)
% 
% Output Variables:
% A - 1 x N vector of fractional area of patches 
% N - number of patches (scalar double)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, N] = DeterministicTilePatches(f, A_tot, A_min)
    
    A_patch = 0;
    A = [];
    N = 0;
    while A_patch < f*A_tot
        if f*A_tot - A_patch >= A_min
            tmp = (f*A_tot - A_patch - A_min)*rand + A_min;
            A = [A tmp/A_tot];
            A_patch = A_patch + tmp;
            N = N + 1;
        else
            break;
        end
    end
    
    if sum(A) > 1
        A = (A_min/A_tot) + sum(A) - 1;
        N = 1;
    end
end