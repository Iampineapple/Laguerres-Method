function [y, q] = Horner(p, r)

%
%************************************************************
%*                                                          *
%*  Routine implementing Horner's scheme for evaluating     *
%*  a polynomial p(x) at a given value r.                   *
%*                                                          *
%*  The degree of p(x) is n and its coefficients are        *
%*  stored in the (n+1)-dimensional vector p                *
%*  (could be column or row):                               *
%*   p(x) = p(1)*x^n + p(2)*x^(n-1) + ... + p(n)*x + p(n+1) *
%*                                                          *
%*                                                          *
%*  Inputs:  p ... coefficient vector (high to low)         *
%*           r ... value to evaluate the polynomial at      *
%*                                                          *
%*  Output:  y ... value of the polynomial at r, p(r)       *
%*           q ... coefficient vector of the polynomial q   *
%*                 such that p(x) = (x-r)*q(x) + p(r)       *
%*                                                          *
%*       Written by Pavel Belik, January 2009               *
%************************************************************
%

n = length(p)-1;   % degree of the polynomial

for i = 2:(n+1)
  p(i) = p(i) + r*p(i-1);   % Overwrite p with q and y
end

q = p(1:n);   % Copy first n entries of p into q
y = p(n+1);   % Last entry of p is y
