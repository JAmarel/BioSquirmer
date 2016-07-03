function mu = HPW_mobility(x)    


Euler = - psi(1);
b1 = 2.74819;
b2 = 0.61465;
c1 = 0.73761;
c2 = 0.52119;



mu = (log(2./x) - Euler + 4*x/pi - (x.^2/2).*log(2./x))./...
                     (1 - x.^3/pi.*log(2./x) + c1*x.^b1./(1+c2*x.^b2) );
                 
                 
end