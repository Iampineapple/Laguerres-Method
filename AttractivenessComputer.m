function [attra] = AttractivenessComputer(x, n, p, dp, ddp, sqrt_poly)
% This function computes the attractiveness of a cycle under Laguerre's Method
%
% This function is called by [attra] = AttractivenessComputer(x, n, p, dp, ddp, sqrt_poly)
%
% Inputs:	x ... a point in the cycle in question
%			n ... the degree of the polynomial, p
%			p ... coefficient vector of the polynomial (high to low)  
%       	dp ... first derivative of p, as coeff. vector 
%         	ddp ... second deriv. of p, as coeff. vector   
%         	sqrt_poly ... expression under the sqrt in     
%               Laguerre's Method 
%
%
% Outputs:	attra ... the attractiveness of the cycle in question
%
% This program was written by Cory Haight-Nali on or around 26 July 2015

%Correct sqrt_poly : 
%sqrt_poly = (n-1)^2*conv(dp, dp) - n*(n-1)*conv(p, ddp);
%sqrt_poly = sqrt_poly(3:end);  % Drop the first two terms, as they will always be zero


%Some initialization
theta = 0;
r = 10^(-4);
h = r * e^(i * theta);


%Calculate L[z] so we can later approximate L'[L[z]]
Lofz = Laguerre_step(x, n, p, dp, ddp, sqrt_poly);

%Approximate L'[z]
deriv = (Laguerre_step(x + h, n, p, dp, ddp, sqrt_poly) - Laguerre_step(x - h, n, p, dp, ddp, sqrt_poly))/(2*h);


%Approximate L'[L[z]]
derivl = (Laguerre_step(Lofz + h, n, p, dp, ddp, sqrt_poly) - Laguerre_step(Lofz - h, n, p, dp, ddp, sqrt_poly))/(2*h);

%Calculate the attraction from the approximations
attra = abs(deriv * derivl);
