function  [out] = Laguerre_step(x, n, p, dp, ddp, sqrt_poly)

%
%************************************************************
%*                                                          *
%*  Routine implementing one step of Laguerre's method      *
%*  for finding a root of a polynomial p(x),                *
%*  given an initial approximation x0.                      *
%*                                                          *
%*  The degree of p(x) is n and its coefficients are        *
%*  passed in the (n+1)-dimensional vector p                *
%*                                                          *
%*                                                          *
%*                                                          *
%*                                                          *
%*  Inputs:  x ... initial approximation to the root        *
%*			 n ... degree of polynomial p                   *
%*           p ... coefficient vector (high to low)         *
%*           dp ... first derivative of p, as coeff. vector *
%*           ddp ... second deriv. of p, as coeff. vector   *
%*           sqrt_poly ... expression under the sqrt in     *
%*               Laguerre's Method                          *   
%*                                                          *
%*  Output:  out ... result of one step of Laguerre's Method*
%*                                                          *
%*       Written by Cory Haight-Nali, May 2015              *
%************************************************************
%


%Correct sqrt_poly : 
%sqrt_poly = (n-1)^2*conv(dp, dp) - n*(n-1)*conv(p, ddp);
%sqrt_poly = sqrt_poly(3:end);  % Drop the first two terms, as they will always be zero

%Calculate value of polynomials at x
px = Horner(p, x);
dpx = Horner(dp, x);
ddpx = Horner(ddp, x);

sqrt_poly_x = Horner(sqrt_poly, x);

%  d = px/dpx * n / (1 + sqrt(sqrt_poly_x/dpx/dpx));

den1 = dpx + sqrt(sqrt_poly_x); %denominator where the plusminus is a plus
den2 = dpx - sqrt(sqrt_poly_x); %denominator where the plusminus is a minus

if (abs(den1) > abs(den2) && abs(den1) > eps) %pick which sign of the plusminus maximizes denominator
	d = n*px / den1; %d is the differential term in Laguerre's method, so z_{k+1} = z_{k} \plusminus d
elseif (abs(den2) > eps)
	d = n*px / den2;
else
	out = x;
	return;
end

out = x - d;

endfunction
