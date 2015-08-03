function [rootsgrid cyclesgridtwo cyclesgridthree cyclesgridfour oldCyclesTwo oldCyclesThree oldCyclesFour oldCyclesTwoMap oldCyclesThreeMap oldCyclesFourMap numroots thecolormap] = BasinLaguerreWithCycles(p, iter, xmin, xmax, ymin, ymax, nx, ny)
%This function is called by BasinLaguerreWithCycles(p, iter, xmin, xmax, ymin, ymax, nx, ny)
%
%This function will take in a polynomial p in coefficant form, high to low, 
%and call Laguerre's method on every point on the rectangle from (xmin, ymin) 
%to (xmax, ymax) with a x resolution of nx pixels, and a y resolution of ny.
%It will then return a graph of the end result of iterating Laguerre's method with up to 
%iter iterations on each of these points.  Each point will be coloured by what root, 
%or cycle it eventually convergea to.
%
%Written by Cory Haight-Nali on or around June 2015

%tol1 is the tolerance for which we match two iterations of Laguerre's method to see if they're the same point
%tol2 is used when comparing a new cycle to old cycles to decide if they are the same cycle
tol1 = 10^(-8);
tol2 = 10^(-4);

% Degree of the polynomial
n = length(p) - 1;

% Roots of the polynomial using the built-in function
rts = roots(p);

%Initialize the initial guesses for the method
z0 = zeros(nx + 1, ny + 1);

re = (0:nx)/nx * (xmax - xmin) + xmin; %set an array for the xcoordinates each pixel (column) starts at
im = (0:ny)/ny * (ymax - ymin) + ymin; %set an array for the ycoordinates each pixel (row) starts at

[xx yy] = meshgrid(re, im);
z0 = (xx + I*yy); %Create a grid of the coordinate points, and set each point in the grid to its coordinates in the complex plane

%Set empty grids to hold what the point converges to
rootsgrid = zeros(ny + 1, nx + 1);
cyclesgridtwo = zeros(ny + 1, nx + 1);
cyclesgridthree = zeros(ny + 1, nx + 1);
cyclesgridfour = zeros(ny + 1, nx + 1);

%Create empty matrices to hold the cycles we find
oldCyclesTwo = [];
oldCyclesThree = [];
oldCyclesFour = [];


%Create empty arrays to hold the mapping of the cycles when we sort them
oldCyclesTwoMap = [];
oldCyclesThreeMap = [];
oldCyclesFourMap = [];



%Start computations
dp = Prime(p);
ddp = Prime(dp);
sqrt_poly = (n-1)^2*conv(dp, dp) - n*(n-1)*conv(p, ddp);
sqrt_poly = sqrt_poly(3:end);  % Drop the first two terms, as they will always be zero


for i = 1:(ny + 1)
	for j = 1:(nx + 1)
		%loop through every point in the plot
		flag = 0;
		foo = 0;
		%and, on each point, call Laguerre_checker on it to see if it has converged to a root or cycle
		%Repeat for iter iterations of Laguerre's method
		while(flag == 0 && foo < iter)
			[one_step two_steps three_steps four_steps flag] = Laguerre_checker(z0(i,j), n, p, dp, ddp, sqrt_poly, tol1);
			z0(i,j) = four_steps;
			foo = foo + 4;
		end
		
		%If no root was found, flag==0, number rootsgrid accordingly
		if (flag == 0)
			rootsgrid(i,j) = n + 2;			
		%If we found a root, flag==1.  Compare the root to the closest precomputed by 
		%rts = roots(p); function
		elseif(flag == 1)
			[aux, index] = min(abs(rts - four_steps));
			%If the iteration stopped moving near a root, and it agrees with 
			%the roots function, mark it as a root according to which root 
			%of the roots function it agrees with.  Otherwise, leave it white.
			if(aux<10^(-4))
				rootsgrid(i,j) = index;
			else
				%Note : this case should never happen
				rootsgrid(i,j) = n + 2;
			end
		%If we have converged to a 2cycle, flag==2.  Get the entire two-cycle, and 
		%see if it's a two-cycle we already have stored.
		%
		%If we already have this two-cycle stored somewhere, mark down the point as 
		%converging to that cycle.  Otherwise, add it to the list of twocycles.
		elseif(flag == 2)
			[cycleFlag index] = LinearCycleSearch([one_step two_steps], oldCyclesTwo, tol2);
			
			%If the cycle is not already saved, calculate the attractiveness,
			%save the cycle and attractiveness in oldCyclesTwo, update oldCyclesTwoMap
			%sort oldCyclesTwo according to attractiveness, 
			%and mark the point in cyclesgridtwo
			if(cycleFlag == 0)
				attra = AttractivenessComputer(one_step, n, p, dp, ddp, sqrt_poly);
				oldCyclesTwo = [oldCyclesTwo; [one_step two_steps attra]];
				oldCyclesTwoMap(end+1) = length(oldCyclesTwoMap) + 1;
				[oldCyclesTwo oldCyclesTwoMap] = CycleSorter(oldCyclesTwo, oldCyclesTwoMap);				
				cyclesgridtwo(i, j) = oldCyclesTwoMap(rows(oldCyclesTwo));
			%else the cycle is an old one we've seen before, mark the point in cyclegridtwo
			else
				cyclesgridtwo(i, j) = oldCyclesTwoMap(index); 		
			end
			rootsgrid(i, j) = n + 2;
		%If we have converged to a three cycle, flag==3.  Get the entire three cycle, and see if we've already stored it.
		%
		%If we already have this 3cycle stored somewhere, mark down the point as converging to that cycle.  Else, add it to the list of 3cycles.
		elseif(flag == 3)
			[cycleFlag index] = LinearCycleSearch([one_step two_steps three_steps], oldCyclesThree, tol2);
			%If the cycle is not already saved, calculate the attractiveness,
			%save the cycle and attractiveness in oldCyclesThree, update oldCyclesThreeMap
			%sort oldCyclesThree according to attractiveness, 
			%and mark the point in cyclesgridthree
			if(cycleFlag == 0)
				attra = AttractivenessComputer(one_step, n, p, dp, ddp, sqrt_poly);
				oldCyclesThree = [oldCyclesThree; [one_step two_steps three_steps attra]];
				oldCyclesThreeMap(end+1) = length(oldCyclesThreeMap) + 1;
				[oldCyclesThree oldCyclesThreeMap] = CycleSorter(oldCyclesThree, oldCyclesThreeMap);				
				cyclesgridthree(i, j) = oldCyclesThreeMap(rows(oldCyclesThree));
			%else the cycle is an old one we've seen before, mark the point in cyclegridthree
			else
				cyclesgridthree(i, j) = oldCyclesThreeMap(index); 		
			end
			rootsgrid(i, j) = n + 2;
			
		%If we have converged to a four cycle, flag==4.  Get the entire four cycle, and see if we've already stored it.
		%
		%If we already have this 4cycle stored somewhere, mark down the point as converging to that cycle.  Else, add it to the list of 4cycles.
		elseif(flag == 4)
			[cycleFlag index] = LinearCycleSearch([one_step two_steps three_steps four_steps], oldCyclesFour, tol2);
			%If the cycle is not already saved, calculate the attractiveness,
			%save the cycle and attractiveness in oldCyclesFour, update oldCyclesFourMap
			%sort oldCyclesFour according to attractiveness, 
			%and mark the point in cyclesgridfour
			if(cycleFlag == 0)
				attra = AttractivenessComputer(one_step, n, p, dp, ddp, sqrt_poly);
				oldCyclesFour = [oldCyclesFour; [one_step two_steps three_steps attra]];
				oldCyclesFourMap(end+1) = length(oldCyclesFourMap) + 1;
				[oldCyclesFour oldCyclesFourMap] = CycleSorter(oldCyclesFour, oldCyclesFourMap);				
				cyclesgridfour(i, j) = oldCyclesFourMap(rows(oldCyclesFour));
			%else the cycle is an old one we've seen before, mark the point in cyclegridfour
			else
				cyclesgridfour(i, j) = oldCyclesFourMap(index); 		
			end
			rootsgrid(i, j) = n + 2;
		else
			rootsgrid(i, j) = n + 2;
		end
	end
end


%Use the maps to create an inverse map to correlate the colors
for i = 1:length(oldCyclesTwoMap)
	inverseOldCyclesTwoMap(oldCyclesTwoMap(i)) = i;
endfor

for i = 1:length(oldCyclesThreeMap)
	inverseOldCyclesThreeMap(oldCyclesThreeMap(i)) = i;
endfor

for i = 1:length(oldCyclesFourMap)
	inverseOldCyclesFourMap(oldCyclesFourMap(i)) = i;
endfor



%Correlate the colors of the two-cycles backwards, according to inversemap
for i = 1:(ny + 1)
	for j = 1:(nx +1)
		if(cyclesgridtwo(i, j) != 0)
			cyclesgridtwo(i, j) = inverseOldCyclesTwoMap(cyclesgridtwo(i,j)) + 1;
		endif
	endfor
endfor


%Correlate the colors of the three-cycles backwards, according to inversemap
for i = 1:(ny + 1)
	for j = 1:(nx +1)
		if(cyclesgridthree(i, j) != 0)
			cyclesgridthree(i, j) = inverseOldCyclesThreeMap(cyclesgridthree(i,j)) + 1;
		endif
	endfor
endfor


%Correlate the colors of the four-cycles backwards, according to inversemap
for i = 1:(ny + 1)
	for j = 1:(nx +1)
		if(cyclesgridfour(i, j) != 0)
			cyclesgridfour(i, j) = inverseOldCyclesFourMap(cyclesgridfour(i, j)) + 1;
		endif
	endfor
endfor


%Shift over the two-cycles onto the rootsgrid, as we'll be plotting the rootsgrid.
num = n + 2; %Number of colors used for the roots
for i = 1:(ny + 1)
	for j = 1:(nx + 1)
		if(cyclesgridtwo(i, j) != 0)
			rootsgrid(i, j) = cyclesgridtwo(i, j) + num;
		end
	end
end

%Shift over the three-cycles onto the rootsgrid, as we'll be plotting the rootsgrid.
num = n + rows(oldCyclesTwo) + 4; %Number of colors used for the roots + number of colors used for the two-cycles.
for i = 1:(ny + 1)
	for j = 1:(nx + 1)
		if(cyclesgridthree(i, j) != 0)
			rootsgrid(i, j) = cyclesgridthree(i, j) + num;
		end
	end
end

%Shift over the four-cycles onto the rootsgrid, as we'll be plotting the rootsgrid.
num = n + rows(oldCyclesTwo) + rows(oldCyclesThree) + 6; %Number of colors used for the roots plus number of colors used for the two-cycles plus number of colors used for the three-cycles.
for i = 1:(ny + 1)
	for j = 1:(nx + 1)
		if(cyclesgridfour(i, j) != 0)
			rootsgrid(i, j) = cyclesgridfour(i, j) + num;
		end
	end
end



%Print the graph
%First, create a matrix using all the appropriate colormaps, so each root/cycle gets mapped with the appropriate color
thecolormap = [hot(n + 2); ocean(int32(rows(oldCyclesTwo)) + 2); gray(int32(rows(oldCyclesThree)) + 2); copper(int32(rows(oldCyclesFour)) + 2)]; 
colormap(thecolormap); 

%Create the graph, and hold it steady to plot the roots and cycles
image([xmin xmax], [ymin ymax], rootsgrid) 
hold on

%Print the roots on the graph as black circles
for i = 1:length(rts)
	plot(real(rts(i)), imag(rts(i)), "ok", "markersize", 4, "linestyle", "none", "markerfacecolor", "black")
end

%Print the 2 cycles on the graph as black + marks
for i = 1:rows(oldCyclesTwo)
	for j = 1:columns(oldCyclesTwo)
		plot(real(oldCyclesTwo(i, j)), imag(oldCyclesTwo(i, j)), "+k", "markersize", 4, "linestyle", "none", "markerfacecolor", "black")
	end
end

%Print the 3 cycles on the graph as black triangles
for i = 1:rows(oldCyclesThree)
	for j = 1:columns(oldCyclesThree)
		plot(real(oldCyclesThree(i, j)), imag(oldCyclesThree(i, j)), "^k", "markersize", 4, "linestyle", "none", "markerfacecolor", "black")
	end
end

%Print the 4 cycles on the graph as black squares
for i = 1:rows(oldCyclesFour)
	for j = 1:columns(oldCyclesFour)
		plot(real(oldCyclesFour(i, j)), imag(oldCyclesFour(i, j)), "sk", "markersize", 4, "linestyle", "none", "markerfacecolor", "black")
	end
end

hold off

%Make sure the values are aligned properly, and make the axis scale equal
axis xy 
axis equal


numroots = n;
endfunction
