function f = mob_perp_approx(x) 

%%% This routine interpolates \alpha_{\perp} response function given in A.Levine and F.McKintosh paper. 
%%% The original function is given in terms
%%% of Struve and Bessel functions. MatLab has serious issues with these
%%% functions: Struve functions are simply not installed, Bessel functions
%%% are unstable for large variables


  
  B1 = 4.1647235912273;
  B2 =  -0.55614445389823;
  C1 = 0.561895056224484;
  C2 = 2.54791114259626;
  D1 = 1.07719599699314;
  D2 = 0.58138830216382;
  D3 = 2.79047959529572;
  D4 = 0.1128659701643252;
  D5 = 1.133498055075156;
  
  f =  (log(2./x) - D1 + 4*x/pi - (x.^2/2).*log(2./x) + D4.*x.^D5 )./...
                     (1 - x.^3/pi.*log(2./x) + C1*x.^B1./(1 + C2*x.^B2) + D2*x.^D3);
  
  
end