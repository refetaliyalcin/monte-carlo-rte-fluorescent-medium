function [sonuc] = k_yagce(x)
%https://www.sciencedirect.com/science/article/abs/pii/S0925838807009577
% x=x*10^9; %m to nm
data_eV=[1  0
    2.3    0
    2.34600760456274e+000	4.92610837438424e-001
2.36882129277567e+000	9.85221674876847e-001
2.39733840304183e+000	1.47783251231527e+000
2.41444866920152e+000	2.46305418719212e+000
2.44296577946768e+000	3.44827586206897e+000
2.46007604562738e+000	4.92610837438424e+000
2.47718631178707e+000	6.40394088669951e+000
2.49429657794677e+000	8.86699507389163e+000
2.50570342205323e+000	1.08374384236453e+001
2.51711026615970e+000	1.37931034482759e+001
2.52851711026616e+000	1.62561576354680e+001
2.54562737642586e+000	2.11822660098522e+001
2.55133079847909e+000	2.36453201970443e+001
2.56273764258555e+000	2.70935960591133e+001
2.57414448669202e+000	3.15270935960591e+001
2.57984790874525e+000	3.64532019704434e+001
2.58555133079848e+000	4.18719211822660e+001
2.59125475285171e+000	4.43349753694581e+001
2.59695817490494e+000	4.67980295566502e+001
2.60266159695817e+000	4.97536945812808e+001
2.60836501901141e+000	5.17241379310345e+001
2.61406844106464e+000	5.36945812807882e+001
2.61977186311787e+000	5.71428571428571e+001
2.62547528517110e+000	5.96059113300493e+001
2.63117870722433e+000	6.20689655172414e+001
2.63688212927757e+000	6.40394088669951e+001
2.64258555133080e+000	6.69950738916256e+001
2.65399239543726e+000	7.14285714285714e+001
2.65969581749049e+000	7.43842364532020e+001
2.66539923954373e+000	7.58620689655172e+001
2.67680608365019e+000	7.73399014778325e+001
2.69391634980989e+000	8.02955665024631e+001
2.69961977186312e+000	8.07881773399015e+001
2.71673003802281e+000	7.98029556650246e+001
2.73384030418251e+000	7.73399014778325e+001
2.75095057034221e+000	7.29064039408867e+001
2.76806083650190e+000	6.69950738916256e+001
2.77946768060837e+000	6.40394088669951e+001
2.78517110266160e+000	5.96059113300493e+001
2.80228136882129e+000	5.46798029556650e+001
2.81368821292776e+000	4.87684729064039e+001
2.83079847908745e+000	4.08866995073892e+001
2.84790874524715e+000	3.30049261083744e+001
2.86501901140684e+000	2.61083743842365e+001
2.87642585551331e+000	2.26600985221675e+001
2.89353612167300e+000	1.77339901477833e+001
2.90494296577947e+000	1.52709359605911e+001
2.92205323193916e+000	1.23152709359606e+001
2.93346007604563e+000	1.08374384236453e+001
2.95627376425856e+000	7.88177339901478e+000
2.97908745247148e+000	5.41871921182266e+000
3.01330798479087e+000	3.44827586206897e+000
3.04752851711027e+000	1.97044334975369e+000
3.09315589353612e+000	1.97044334975369e+000
3.18441064638783e+000	1.97044334975369e+000
3.28136882129278e+000	1.97044334975369e+000
3.38403041825095e+000	1.97044334975369e+000
3.43536121673004e+000	2.95566502463054e+000
3.44676806083650e+000	3.44827586206897e+000
3.46387832699620e+000	4.92610837438424e+000
3.48098859315589e+000	6.40394088669951e+000
3.49809885931559e+000	7.88177339901478e+000
3.50950570342205e+000	1.03448275862069e+001
3.51520912547529e+000	1.13300492610837e+001
3.53231939163498e+000	1.42857142857143e+001
3.54372623574144e+000	1.67487684729064e+001
3.55513307984791e+000	1.92118226600985e+001
3.57224334600760e+000	2.26600985221675e+001
3.57794676806084e+000	2.41379310344828e+001
3.58935361216730e+000	2.66009852216749e+001
3.60076045627376e+000	2.80788177339901e+001
3.61216730038023e+000	2.90640394088670e+001
3.62927756653992e+000	3.00492610837438e+001
3.64068441064639e+000	3.00492610837438e+001
3.65209125475285e+000	2.95566502463054e+001
3.67490494296578e+000	2.80788177339901e+001
3.68631178707224e+000	2.70935960591133e+001
3.69771863117871e+000	2.56157635467980e+001
3.72053231939163e+000	2.26600985221675e+001
3.74334600760456e+000	1.92118226600985e+001
3.76045627376426e+000	1.62561576354680e+001
3.77186311787072e+000	1.47783251231527e+001
3.79467680608365e+000	1.23152709359606e+001
3.82319391634981e+000	9.85221674876847e+000
3.84600760456274e+000	8.37438423645320e+000
3.88022813688213e+000	6.89655172413793e+000
3.90304182509506e+000	6.40394088669951e+000
3.97148288973384e+000	6.40394088669951e+000
4.03992395437262e+000	7.38916256157636e+000
4.08555133079848e+000	7.38916256157636e+000
4.13117870722433e+000	8.37438423645320e+000
4.15969581749049e+000	8.37438423645320e+000
4.19961977186312e+000	9.85221674876847e+000];
% data_eV(:,1)=data_eV(:,1);
data_eV(:,2)=data_eV(:,2)*20*1000/max(data_eV(:,2));
data_nm=data_eV;
data_nm(:,1)=12.4 * 10^-7 ./data_eV(:,1);
data_nm=flip(data_nm);
data_nm(:,2)=data_nm(:,1).*data_nm(:,2)/(4*pi);

sonuc=interp1(data_nm(:,1),data_nm(:,2),x);
end

