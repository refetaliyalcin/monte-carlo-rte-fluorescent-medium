notOctave = exist('OCTAVE_VERSION', 'builtin') == 0;


beta_non_fl=zeros(number_wl,1);
beta_fl=zeros(number_wl,1);
%abs and scat coef


r_vector=(radius-2.45*sigma):sigma/200:(radius+2.45*sigma);
weight_vector = normpdf(r_vector,radius,sigma);
weight_vector=weight_vector/trapz(r_vector,weight_vector);


mu_tot_arr_sigma=zeros(length(lamda),length(r_vector));
beta_sigma=zeros(length(lamda),length(r_vector));
alfa_sigma=zeros(length(lamda),length(r_vector));
scat_prob_arr_sigma=zeros(length(lamda),length(r_vector));
g_arr_sigma=zeros(length(lamda),length(r_vector));
QY_modified_sigma=zeros(length(lamda),length(r_vector));
for z=1:length(r_vector)
    Area=pi*r_vector(z)^2;
    V=4*pi*r_vector(z)^3/3;
    for i=1:number_wl
        x=2*pi*r_vector(z)*n_medium(i)/lamda(i);
        m=(n_phosphor(i)+1i*k_phosphor(i))/n_medium(i);
        fonksiyon=Mie(m,x);
        Qabs=fonksiyon(3);
        Qsca=fonksiyon(2);
        g_arr_sigma(i,z)=fonksiyon(5);
        Qsca=Qsca*(1.06-g_arr_sigma(i,z));
        Qabs=Qabs*1.47;
        alfa_sigma(i,z)=f_v*Qsca*Area/V;
        beta_fl=f_v*Qabs*Area/V;
        beta_non_fl=(1-f_v)*4*pi*k_medium(i)/lamda(i); %absorption by medium (non fluorescent)
        beta_sigma(i,z)=beta_fl+beta_non_fl;
        QY_modified_sigma(i,z)=QY(i)*beta_fl/beta_sigma(i,z); %Modified version of QY is used since the probability of re-emitting the absorbed rays by non fluorescent part is zero.
        mu_tot_arr_sigma(i,z)=alfa_sigma(i,z)+beta_sigma(i,z);
        scat_prob_arr_sigma(i,z)=alfa_sigma(i,z)/mu_tot_arr_sigma(i,z);
    end
end
ext_tot=trapz(r_vector,(weight_vector.*mu_tot_arr_sigma)');
scat_prob=trapz(r_vector,(weight_vector.*scat_prob_arr_sigma)');
QY_modified=trapz(r_vector,(weight_vector.*QY_modified_sigma)');
g=trapz(r_vector,(weight_vector.*g_arr_sigma.*scat_prob_arr_sigma)')./trapz(r_vector,(weight_vector.*scat_prob_arr_sigma)');
beta=trapz(r_vector,(weight_vector.*beta_sigma)');
alfa=trapz(r_vector,(weight_vector.*alfa_sigma)');
% pdf ve cdf

data_flo=flo_emission_data();
flo_wl_start=data_flo(1,1);
flo_wl_end=data_flo(end,1);

data_1=data_flo(:,1);
data_2=data_flo(:,2);
wl_arr=(flo_wl_start:flo_wl_end)';
result_arr=interp1(data_1,data_2,wl_arr);
wave_flo=flo_wl_start:flo_wl_end;
wave_flo_no=length(wave_flo);
x=zeros(wave_flo_no,1);
y=zeros(wave_flo_no,1);
cdf=zeros(wave_flo_no,1);

for i=1:wave_flo_no
   x(i)=(i-1)/(wave_flo_no-1);
   y(i)=result_arr(i);
end
y=y/trapz(x,y);
for i=2:wave_flo_no
    cdf(i)=trapz(x(1:i),y(1:i));
end

cdf_fl_new=zeros(10001,1);
inv_cdf=zeros(10001,1);
for i=1:10001
   cdf_fl_new(i)=i-1;
   inv_cdf(i)=interp1(cdf,wave_flo,cdf_fl_new(i)/10000);
end


fig1=figure(1);

plot(wl,exc_yagce(wl),'-k',data_1,data_2/max(data_2),'--k','LineWidth',2)
ylabel('Normalized Intensity [a.u.]')
xlabel('Wavelength [nm]')
xlim([400 750])
legend('Excitation','Emission','Location','NorthEast')
saveas(fig1,'abs_emis.fig')
saveas(fig1,'abs_emis.emf')

% refractive indices
if notOctave %yyaxis problem in octave
    fig2=figure(2);
    set(fig2,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
    hold on
    yyaxis left
    plot(wl,n_phosphor,'-k','LineWidth',2)
    ylabel('n')
    ylim([1.7 2.1])
    yyaxis right
    plot(wl,k_phosphor,'--k','LineWidth',2)
    ylabel('k')
    hold off
    box on
    xlabel('Wavelength [nm]')
    xlim([400 end_wl])
    legend('n','k','Location','SouthEast')
    saveas(fig2,'ref_ind.fig')
    saveas(fig2,'ref_ind.emf')

%     absorption and scat coeff
    fig3=figure(3);
    set(fig3,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
    hold on
    yyaxis left
    plot(wl,alfa,'-k',wl,beta,':k','LineWidth',2)
    ylabel('Scattering and Absorption Coefficient [1/m]')
%     ylim([0  1.1*max(alfa)])
    yyaxis right
    plot(wl,g,'--k','LineWidth',2)
    ylabel('Asymmetry Parameter (g)')
    ylim([0.85 1])
    hold off
    box on
    xlabel('Wavelength [nm]')
    xlim([400 end_wl])
    legend('Scattering Coefficient','Absorption Coefficient','Asymmetry Parameter','Location','SouthEast')
    saveas(fig3,'scat_abs_coef.fig')
    saveas(fig3,'scat_abs_coef.emf')

%     CDF and PDF
    fig5=figure(5);
    set(fig5,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
    hold on
    yyaxis left
    plot(flo_wl_start:flo_wl_end,y/trapz(flo_wl_start:flo_wl_end,y),'-k','LineWidth',2);
    xlabel('Wavelength [nm]')
    ylabel('PDF [1/nm]')
    ylim([0 0.01])
    xlim([flo_wl_start flo_wl_end])
    yyaxis right
    plot(flo_wl_start:flo_wl_end,cdf,'--k','LineWidth',2);
    xlabel('Wavelength [nm]')
    ylabel('CDF')
    ylim([0 1])
    xlim([flo_wl_start flo_wl_end])
    legend('Probability Density Function','Cumulative Distribution Function','Location','northwest')
    hold off
    box on
    saveas(fig5,'pdf_cdf.fig')
    saveas(fig5,'pdf_cdf.emf')
end
%inv cdf
fig6=figure(6);
plot(linspace(0,1,length(inv_cdf)),inv_cdf,'-k','LineWidth',2);
ylabel('Wavelength [nm]')    
xlabel('Random Number')    
xlim([0 1])
ylim([flo_wl_start flo_wl_end])
legend('Inverse of Cumulative Distribution Function','Location','southeast')
saveas(fig6,'invcdf.fig')
saveas(fig6,'invcdf.emf')
    