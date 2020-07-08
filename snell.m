function alive=snell(s_z,n_medium,k_medium,n_subs,k_subs)
% n_medium,k_medium,n_subs,k_subs
    cos_teta=abs(s_z);
    if (s_z>0) 
       n_outside=n_subs; 
       k_outside=k_subs; 
    else
       n_outside=1; %air
       k_outside=0;
    end
    if (n_outside==n_medium && k_outside==k_medium) % n'ler ayniysa pasaport sorma dogrudan geç
        reflectance=0;
    elseif (cos_teta>0.9999) %dik gelmissen açiyi karistirma
        reflectance=((n_medium-n_outside)^2+(k_medium-k_outside)^2)/((n_medium+n_outside)^2+(k_medium+k_outside)^2);
    elseif (cos_teta<(0.0001)) %paralel gelmissen yansit
        reflectance=1;
    else
        % genel durum bu howell 5th ed. sahife 741
        sin_teta=sqrt(1-cos_teta*cos_teta);
        carpan2=(n_medium-1i*k_medium)/(n_outside-1i*k_outside);
        sin_x2=sin_teta*carpan2;
        cos_x2=sqrt(1-sin_x2*sin_x2);
        E_parallel=(cos_teta/cos_x2-carpan2)/(cos_teta/cos_x2+carpan2);
        R_parallel=E_parallel*conj(E_parallel);
        E_orth=-(cos_x2/cos_teta-carpan2)/(cos_x2/cos_teta+carpan2);
        R_orth=E_orth*conj(E_orth);
        reflectance=real(R_parallel+R_orth)*0.5;
    end
    Rastgele=rand();
    
    if Rastgele>reflectance
        %pass
        alive=0;
    else
        %reflect
        alive=1;
    end
 end 