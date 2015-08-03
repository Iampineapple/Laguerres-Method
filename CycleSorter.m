function [cycles map] = CycleSorter(cycles, map)
%This function sorts a matrix of cycles where the final column of the matrix
%is the attractiveness of the cycles.  It also correspondingly re-arranges the map,
%so we know which cycle to color which.
%
% This function is called by [cycles map] = CycleSorter(cycles, map)
%
% Inputs:	cycles ... a matrix of cycles, augmented by the attractiveness of the cycle
%					the final cycle needs to be sorted into the matrix
%			map  ... array of the color corresponding to the cycle in cycles(row)
%
% Outputs:	cycles ... a matrix of cycles, augmented by the attractiveness of the cycle
%			map  ... array of the color corresponding to the cycle in cycles(row)
%
% This program was written by Cory Haight-Nali on or around 27 July 2015

j = rows(cycles);


if(j>=2)
	while(j>1 && cycles(j-1, columns(cycles)) > cycles(j, columns(cycles)))
		%Swap the cycles rows
		temp = cycles(j-1, 1:columns(cycles));
		cycles(j-1, 1:columns(cycles)) = cycles(j, 1:columns(cycles));
		cycles(j, 1:columns(cycles)) = temp;
		
		%Swap the corresponding map elements
		temp = map(j-1);
		map(j-1) = map(j);
		map(j) = temp;
		
		%decrement j
		j = j -1;
	endwhile
endif


	