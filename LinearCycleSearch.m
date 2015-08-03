function [flag index] = LinearCycleSearch(newCycle, oldCycles, tol)
%called by [flag index] = LinearCycleSearch(newCycle, oldCycles, tol)
%
% newCycle is an array of all the elements of the new cycle 
% of which we will check inclusion in oldCycle.  oldCycles is a 2d array consisting 
% of arrays of old, previously dealt with cycles
%
% flag == 0 if the items in newCycle were not found in oldCycles
% flag == 1 if one of the items in newCycle was found in oldCycles
%
% Index indicates the position in oldCycles where the item in newCycle was found;
% that is, index indicates which saved cycle newCycle was found in.  
% If newCycle was not found, index will be return zero.
%
% This program was written by Cory Haight-Nali on or around June 2015

flag = 0;
index = 0;

%if oldCycles is an empty matrix, quit
if(length(oldCycles) == 0)
	return;
end

for newCycleIndex = 1:length(newCycle)
	for oldCycleIndex = 1:rows(oldCycles)
		for oldCyclePart = 1:columns(oldCycles)-1
				if( abs(newCycle(newCycleIndex) - oldCycles(oldCycleIndex, oldCyclePart)) < tol * oldCycles(oldCycleIndex, oldCyclePart))
				flag = 1;
				index = oldCycleIndex;
				return;
			end
		end
	end
end

endfunction