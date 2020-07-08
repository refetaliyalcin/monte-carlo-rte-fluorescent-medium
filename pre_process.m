notOctave = exist('OCTAVE_VERSION', 'builtin') == 0;

alfa=zeros(number_wl,1);
beta_non_fl=zeros(number_wl,1);
beta_fl=zeros(number_wl,1);
beta=zeros(number_wl,1);
g=zeros(number_wl,1);
QY_modified=zeros(number_wl,1);

abs_coeff=abs_coeff_yagce(lamda);

k=2*pi./lamda;
k_yc=abs_coeff./(2*k);
n_yc=n_yagce(lamda);

V=4*pi*radius^3/3;

%abs and scat coef
for i=1:number_wl
    ev_pigment_n=n_yc(i);
    ev_pigment_k=k_yc(i);
    x=2*pi*radius*n_medium(i)/lamda(i);
    m=(ev_pigment_n+1i*ev_pigment_k)/n_medium(i);
    fonksiyon=Mie(m,x);
    Qext=fonksiyon(1);
    Qsca=fonksiyon(2);
    Csca=pi*radius^2*Qsca;
    Cext=pi*radius^2*Qext;
    Cabs=Cext-Csca;
    alfa(i)=f_v*Csca/V;
    beta_fl(i)=f_v*Cabs/V; %absorption by fluorescence
    beta_non_fl(i)=4*pi*k_medium(i)/lamda(i); %absorption by medium (non fluorescent)
    beta(i)=beta_fl(i)+beta_non_fl(i);
    QY_modified(i)=QY*beta_fl(i)/beta(i); %Modified version of QY is used since the probability of re-emitting the absorbed rays by non fluorescent part is zero.
    g(i)=fonksiyon(5);
end

ext_tot=alfa+beta;
scat_prob=alfa./ext_tot;


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

plot(wl,abs_coeff/max(abs_coeff),'-k',data_1,data_2/max(data_2),'--k','LineWidth',2)
ylabel('Normalized Intensity')
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
    plot(wl,n_yc,'-k','LineWidth',2)
    ylabel('n')
    ylim([1.7 2.1])
    yyaxis right
    plot(wl,k_yc,'--k','LineWidth',2)
    ylabel('k')
    hold off
    box on
    xlabel('Wavelength [nm]')
    xlim([400 1000])
    legend('n','k','Location','SouthEast')
    saveas(fig2,'ref_ind.fig')
    saveas(fig2,'ref_ind.emf')

%     absorption and scat coeff
    fig3=figure(3);
    set(fig3,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
    hold on
    yyaxis left
    plot(wl,alfa,'-k',wl,beta*10,'--k','LineWidth',2)
    ylabel('Scattering and Absorption Coefficient [1/m]')
    ylim([0  1.4*max(alfa)])
    yyaxis right
    plot(wl,g,':k','LineWidth',2)
    ylabel('Asymmetry Parameter (g)')
    ylim([0.75 1])
    hold off
    box on
    xlabel('Wavelength [nm]')
    xlim([400 1000])
    legend('Scattering Coefficient','Absorption Coefficient x 10','Asymmetry Parameter','Location','NorthEast')
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
    