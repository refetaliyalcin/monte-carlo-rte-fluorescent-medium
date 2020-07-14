% Codded by Refet Ali YALCIN. You can change and distribute the code
% but keep this commented part. The codes come WITHOUT ANY WARRANTY 
% In case of use, please cite the articles as:
% https://doi.org/10.1088/2053-1591/ab28b8 Improving photosynthetic efficiency using greenhouse coatings with scattering and fluorescent pigments
% https://doi.org/10.1016/j.biosystemseng.2020.02.007 Improving crop production in solar illuminated vertical farms using fluorescence coatings
function [absorption_no,reflect_no,trans_no]=monte_carlo(h,wl,scat_prob,ext_tot,g,QY_modified,start_wl,number_wl,inv_cdf,cos_teta_prime,sur_reflection,n_medium,k_medium,n_subs,k_subs)
    absorption_no=0;
    reflect_no=zeros(number_wl,1);
    trans_no=zeros(number_wl,1);
    wl_index=wl-start_wl+1;
    reflect_no(wl_index)=sur_reflection; 
    energy=1-sur_reflection;
    x=0; %position vector x component
    y=0; %position vector y component
    z=0; %position vector z component
    phi=2*pi*rand;
    sin_teta_prime=sqrt(1-cos_teta_prime*cos_teta_prime);
    sin_phi=sin(phi);
    cos_phi=sqrt(1-sin_teta_prime*sin_teta_prime);
    s_x=sin_teta_prime*sin_phi; %direction vector x component
    s_y=sin_teta_prime*cos_phi; %direction vector y component
    s_z=cos_teta_prime; %direction vector z component
    alive=1; %1 if ray bundle alive, 0 if dead
    l_beta=-log(rand())/ext_tot(wl_index); %extinction length
    while alive
        if (s_z>0) %if bundle is moving in z direction to up, calculate distance to upper boundary
            l_w = (h - z)/s_z; %distance to upper boundary
        else %if bundle is moving in -z direction, calculate distance to lower boundary
            l_w = -z/s_z; %distance to upper boundary
        end
        if l_w<l_beta %if distance to wall is smaller than distance to extinction, bundle hits boundary
            min_index=1;
            min_l=l_w;
        else %if distance to wall is greater than distance to extinction, bundle extincts
            min_index=2;
            min_l=l_beta;
        end
        x=x+min_l*s_x;%move the ray bundle in x
        y=y+min_l*s_y;%move the ray bundle in y
        z=z+min_l*s_z;%move the ray bundle in z
        if (min_index==1) % ray bundle reaches boundary
            alive=snell(s_z,n_medium(wl_index),k_medium(wl_index),n_subs(wl_index),k_subs(wl_index));
            if (alive==0) %ray bundle escapes
                if s_z>0
                    trans_no(wl_index)=trans_no(wl_index)+energy;
                else    
                    reflect_no(wl_index)=reflect_no(wl_index)+energy;
                end
            else %ray bundle reflects
                s_z=-s_z;
                l_beta=l_beta-l_w;
            end
        else % ray bundle extinct
            if rand()<scat_prob(wl_index)  % ray bundle scatters
                [s_x,s_y,s_z]=scatter(g(wl_index),s_x,s_y,s_z);
                l_beta=-log(rand())/ext_tot(wl_index);
            else % ray bundle absorbed
                if rand()>QY_modified(wl_index) %check if ray bundle is absorbed for good
                    alive=0;
                    absorption_no=absorption_no+energy;
                else %ray is reemitted after absorption
                    angles=random_angles(); % get sin_theta,cos_theta,sin_phi,cos_phii for random emission
                    s_x=angles(1);
                    s_y=angles(2);
                    s_z=angles(3);
                    wl_new=round(inv_cdf(ceil(rand()*10000)));% find new wavelength of the ray bundle
                    absorption_no=absorption_no+energy*(1-wl/wl_new); % energy change in host medium due to fluorescent conversion
                    energy=energy*wl/wl_new;  % energy change in ray bundle due to fluorescent conversion
                    %gather and recalculate properties for new lamda
                    wl=wl_new;
                    wl_index=wl-start_wl+1;
                    l_beta=-log(rand())/ext_tot(wl_index);
                end
            end
        end
    end
end