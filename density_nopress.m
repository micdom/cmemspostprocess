% M. De Dominicis (2015)
% density according to the Jackett et al. (2006)
% Algorithms for Density, Potential Temperature, Conservative Temperature, and the Freezing Temperature of Seawater
% th : Potential temperature (degrees Celsius)
% s: Salinity (PSU)
% p: pressure
% Check values
% PEA_density(35,25,2000)=1031.65056056576 
% PEA_density(20,20,1000)=1017.72886801964
% PEA_density(40,12,8000)=1062.95279820631
function [dens] = density_nopress(s,th)
format long
p=0

th2 =th.*th;
sqrts = s.^(1/2);
pth = p.*th;
anum = 9.9984085444849347e+02 + th.*( 7.3471625860981584e+00 + th.*(-5.3211231792841769e-02 + th.*3.6492439109814549e-04)) + ...
            s.*( 2.5880571023991390e+00 - th.*6.7168282786692355e-03 + s.*1.9203202055760151e-03) + ...
            p.*(1.1798263740430364e-02 + th2.*9.8920219266399117e-08 + s.*4.6996642771754730e-06 - p.*(2.5862187075154352e-08 + th2.*3.2921414007960662e-12));

aden = 1.0 + th.*(7.2815210113327091e-03 + th.*(-4.4787265461983921e-05 + th.*(3.3851002965802430e-07 + th.*1.3651202389758572e-10))) + ...
            s.*(1.7632126669040377e-03 - th.*(8.8066583251206474e-06 + th2.*1.8832689434804897e-10) + sqrts.*(5.7463776745432097e-06 + th2.*1.4716275472242334e-09)) + ...    
            p.*(6.7103246285651894e-06 - pth.*(th2*2.4461698007024582e-17 + p.*9.1534417604289062e-18));

dens = anum./aden;
 

      
 

 