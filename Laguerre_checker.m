function [one_step two_steps three_steps four_steps flag] = Laguerre_checker(z, n, p, dp, ddp, sqrt_poly, tol)
%This function checks if L(z) ~= z, L(L(z)) ~= z, and if L^4(z) ~= z
%
%It is called by [one_step two_steps three_steps four_steps flag] 
%		= Laguerre_checker(z, n, p, dp, ddp, sqrt_poly, tol)
%
% Inputs: z ... initial point
%         p ... coefficient vector of the polynomial (high to low)         
%         dp ... first derivative of p, as coeff. vector 
%         ddp ... second deriv. of p, as coeff. vector   
%         sqrt_poly ... expression under the sqrt in     
%               Laguerre's Method   
%		  tol ... minimum tolerance for which x and the iterated 
%				Laguerre's Method can differ by and be considered
%				the same                       
%
% Outputs: one_step ... output point after one step of Laguerre's Method
%		   two_steps ... output point after two steps of Laguerre's Method
%		   three_steps ... output point after three steps of Laguerre's Method
%		   four_steps ... output point after four steps of Laguerre's Method
%          flag ... flag to return what was found
%					flag == 0 means that neither was true
%					flag == 1 means that L(z) ~= z, that we have a root
%					flag == 2 means that L(L(z)) ~= z (and thus also L^(4) ~= z)
%					flag == 3 means that L^3(z) ~= z, so we have a 3 cycle
%					flag == 4 means that L^4(z) ~= z, so we have a 4 cycle
%
%Written by Cory Haight-Nali on or around 20 May 2015

%Correct sqrt_poly : 
%sqrt_poly = (n-1)^2*conv(dp, dp) - n*(n-1)*conv(p, ddp);
%sqrt_poly = sqrt_poly(3:end);  % Drop the first two terms, as they will always be zero


flag = 0;

one_step = Laguerre_step(z, n, p, dp, ddp, sqrt_poly);
two_steps = Laguerre_step(one_step, n, p, dp, ddp, sqrt_poly);
three_steps = Laguerre_step(two_steps, n, p, dp, ddp, sqrt_poly);
four_steps = Laguerre_step(three_steps, n, p, dp, ddp, sqrt_poly);


if(abs(z - one_step) < (z * tol))
	flag = 1;
	return
elseif(abs(z - two_steps) < (z * tol))
	flag = 2;
	return
elseif(abs(z - three_steps) < (z * tol))
	flag = 3;
	return
elseif(abs(z - four_steps) < (z * tol))
	flag = 4;
	return
end

endfunction