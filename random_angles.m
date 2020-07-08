function [ s ] = random_angles()
    costheta = 1.0 - 2.0*rand();
    sintheta = sqrt(1.0 - costheta*costheta);
    psi = 2.0*pi*rand();
    cospsi = cos(psi);
    if (psi < pi)
        sinpsi = sqrt(1.0 - cospsi*cospsi); 
    else
        sinpsi = -sqrt(1.0 - cospsi*cospsi);
    end    
    s(1) = sintheta*cospsi;
    s(2) = sintheta*sinpsi;
    s(3) = costheta;
end 