function result=n_yagce(lam)
% https://doi.org/10.1364/AO.49.000247
lam=lam*10^6; %m to um
result=sqrt(2.08745+1.2081*lam.^2./(lam.^2-0.02119) +17.2049*lam.^2./(lam.^2-1404.45));
