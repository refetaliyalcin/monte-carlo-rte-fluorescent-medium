% Codded by Refet Ali YALCIN. You can change and distribute the code
% please keep this part. The codes come WITHOUT ANY WARRANTY 
% In case of use, following articles that uses the code can be cited:
% https://doi.org/10.1088/2053-1591/ab28b8 Improving photosynthetic efficiency using greenhouse coatings with scattering and fluorescent pigments
% https://doi.org/10.1016/j.biosystemseng.2020.02.007 Improving crop production in solar illuminated vertical farms using fluorescence coatings
clear all
close all
clc

start_wl=400; % starting wavelength in nm, must be an integer
end_wl=1000; % last wavelength of the area of interest, must be an integer

repeat_no=50000; % # of montecarlo simulations for each wavelength
h=1*10^-3; %thickness of coating in meters
radius=5000*10^-9; % radius of fluorescent particles in meters
f_v=0.01;  % volume fraction of phosphor particles
QY=0; %quantum yield

polar_angle=linspace(0,89.99999,30); %polar angle is set but not used. you can loop the code over polar angle if required
polar_angle_rad=polar_angle*pi/180;
wl=(start_wl:end_wl)';
number_wl=length(wl);
lamda=wl*10^-9;

n_medium=PMMA_n(lamda);
k_medium=PMMA_k(lamda);
% k_medium=zeros(number_wl,1); % enable for non absorbing medium case
n_subs=ones(number_wl,1); %substrate is air
k_subs=zeros(number_wl,1);
pre_process % call pre process to calculate coefficients

% Below part calculates reflectance and refraction at the air - medium interface
cos_teta_prime=zeros(length(lamda),length(polar_angle)); %the cos of the ray after refracted from air to medium
sur_reflection=zeros(length(lamda),length(polar_angle)); %surface reflection at air-medium interfece
for o=1:length(polar_angle_rad)
    for i=1:length(lamda)
        cos_teta_prime(i,o)=cos(F_fresnel_2(n_medium(i),k_medium(i),polar_angle_rad(o)));
        cos_teta=cosd(polar_angle(o));
        sin_teta=sqrt(1-cos_teta*cos_teta);
        carpan2=1/(n_medium(i)-1i*k_medium(i));
        sin_x2=sin_teta*carpan2;
        cos_x2=sqrt(1-sin_x2*sin_x2);
        carpan1=cos_teta/cos_x2;
        carpan3=cos_x2/cos_teta;
        E_parallel=(carpan1-carpan2)/(carpan1+carpan2);
        R_parallel=E_parallel*conj(E_parallel);
        E_orth=(carpan3-carpan2)/(carpan3+carpan2);
        R_orth=E_orth*conj(E_orth);
        reflectance=real(R_parallel+R_orth)*0.5;
        sur_reflection(i,o)=reflectance;
    end
end
o=1; %for normal incidance;
% end of refraction and surface reflection calculation

% set database for photons
db_absorption_no=zeros(number_wl,1);
db_reflect_no=zeros(number_wl,number_wl);
db_trans_no=zeros(number_wl,number_wl);
tic

for k=start_wl:end_wl
    absorption_no=0;
    reflect_no=zeros(number_wl,1);
    trans_no=zeros(number_wl,1);
%     surface_r_no=0;
    wl_index=k-start_wl+1;
    for i=1:repeat_no % use parallel calculation
        [absorption_no_new,reflect_no_new,trans_no_new] = monte_carlo(h,k,scat_prob,ext_tot,g,QY_modified,start_wl,number_wl,inv_cdf,cos_teta_prime(wl_index,o),sur_reflection(wl_index,o),n_medium,k_medium,n_subs,k_subs);
        absorption_no=absorption_no + absorption_no_new;
        reflect_no=reflect_no + reflect_no_new;
        trans_no=trans_no + trans_no_new;
    end
    db_reflect_no(:,wl_index)=reflect_no;
    db_absorption_no(wl_index)=absorption_no;
    db_trans_no(:,wl_index)=trans_no;
    
    clc
    disp([num2str(floor(wl_index*100/number_wl)),'% has been completed.']);
end
toc
disp(['100% has been completed.']);
post_process
