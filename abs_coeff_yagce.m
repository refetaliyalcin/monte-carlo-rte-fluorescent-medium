function [sonuc] = abs_coeff_yagce(x)
%data from doi:10.1016/j.jlumin.2006.05.004
x=x*10^9; %m to nm
data=[3.78E+02	0.0000E+00
3.92E+02	0.0000E+00
4.01E+02	1.3667E+00
4.09E+02	6.8337E+00
4.15E+02	1.2756E+01
4.29E+02	8.5649E+01
4.39E+02	1.6036E+02
4.46E+02	1.6264E+02
4.74E+02	1.6264E+02
4.80E+02	1.5809E+02
5.03E+02	2.1868E+01
5.05E+02	1.5945E+01
5.11E+02	6.8337E+00
5.21E+02	1.3667E+00
5.27E+02	0.0000E+00
1000.1 0];
conv_cm1_m1=100;
sonuc=interp1(data(:,1),conv_cm1_m1*data(:,2),x);
end

