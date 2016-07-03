function f = mob_par_approx(x)

%%% This routine interpolates \alpha_{\|} response function given in
%%% A.Levine and F.McKintosh paper. The original function is given in terms
%%% of Struve and Bessel functions. MatLab has serious issues with these
%%% functions: Struve functions are simply not installed, Bessel functions
%%% are unstable for large variables

  B1 = 2.573241064;
  B2 =  6.6283*10^(-6);
  C1 = 23.674674998465;
  C2 = 22.382432191646;
  D1 = 0.07726822808401;
  D2 = 0.267832702242734;
  D3 = 1.190656469633888;
  D4 = 0.67720049611942;
  D5 = 2.02475492047006;
  
  f =  (log(2./x) - D1 + 4*x/pi - (x.^2/2).*log(2./x) + D4.*x.^D5 )./...
                     (1 - x.^3/pi.*log(2./x) + C1*x.^B1./(1 + C2*x.^B2) + D2*x.^D3);
end

 