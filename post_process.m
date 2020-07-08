
BB_Source=zeros(number_wl,1);
for i=1:number_wl
    BB_Source(i)=I_bb(lamda(i),5800);
end

BB_Source=1000*BB_Source/sum(BB_Source); %normalize source at 1000W/m^2

% error('123');
prop_t=db_trans_no/repeat_no;
prop_r=db_reflect_no/repeat_no;
prop_a=db_absorption_no/repeat_no;

total_photon=sum(sum(db_reflect_no,2))+sum(sum(db_trans_no,2))+sum(db_absorption_no);
one=total_photon/(repeat_no*number_wl); %this should be one for check
q_ref_lamda=zeros(number_wl,1);
q_tra_lamda=zeros(number_wl,1);
q_abs_lamda=BB_Source.*prop_a;

for i=start_wl:end_wl
    wl_index=i-start_wl+1;
    q_ref_lamda(wl_index)=sum(BB_Source'.*(prop_r(wl_index,:)));
    q_tra_lamda(wl_index)=sum(BB_Source'.*(prop_t(wl_index,:)));
end
ref_lamda=q_ref_lamda./BB_Source;


%spectral fluxes
fig7 = figure(7);
plot(wl,BB_Source,'--k',wl,q_ref_lamda,':k',wl,q_abs_lamda,'-k',wl,q_tra_lamda,'-.k','LineWidth',2)
xlabel('Wavelength [nm]') % x-axis label
ylabel('Spectral Flux [W/m^2\cdotnm]') % y-axis label
ylim([0 2.5])
xlim([400 1000])
box on
legend('Incoming Flux','Reflected Flux','Absorbed Flux','Transmitted Flux','Location','northeast')
saveas(fig7,'Spectral_Flux.fig')
saveas(fig7,'Spectral_Flux.emf')

if notOctave %yyaxis problem in octave
    %2d transmittance
    fig8=figure(8);
    log10_t=log10(prop_t);
    enmini=log10(min(prop_t(prop_t>0)));
    log10_t(~isfinite(log10_t)) = enmini;
    contourf(wl,wl,log10_t,'edgecolor','none')
    xlabel('Excitation Wavelength [nm]') % x-axis label
    ylabel('Emission Wavelength [nm]') % y-axis label
    xlim([400 550])
    ylim([450 730])
    h=colorbar;
    colormap(flipud(gray))
    ylabel(h, 'log(T_{\lambda,\lambda''})')
    daspect([1 1 1])
    saveas(fig8,'trans_2d.fig')
    saveas(fig8,'trans_2d.emf')

    fig11=figure(11);
    qin_450=BB_Source(find(wl==450));
    qin_500=BB_Source(find(wl==500));
    qin_550=BB_Source(find(wl==550));
    qin_800=BB_Source(find(wl==800));
    qin_all=sum(BB_Source);
    qr_450=q_ref_lamda(find(wl==450));
    qr_500=q_ref_lamda(find(wl==500));
    qr_550=q_ref_lamda(find(wl==550));
    qr_800=q_ref_lamda(find(wl==800));
    qr_all=sum(q_ref_lamda);
    qt_450=q_tra_lamda(find(wl==450));
    qt_500=q_tra_lamda(find(wl==500));
    qt_550=q_tra_lamda(find(wl==550));
    qt_800=q_tra_lamda(find(wl==800));
    qt_all=sum(q_tra_lamda);
    qa_450=q_abs_lamda(find(wl==450));
    qa_500=q_abs_lamda(find(wl==500));
    qa_550=q_abs_lamda(find(wl==550));
    qa_800=q_abs_lamda(find(wl==800));
    qa_all=sum(q_abs_lamda);
    qtot_450=qr_450+qt_450+qa_450;
    qtot_500=qr_500+qt_500+qa_500;
    qtot_550=qr_550+qt_550+qa_550;
    qtot_800=qr_800+qt_800+qa_800;
    qtot_all=qr_all+qt_all+qa_all;

    qName = {'Incoming';'R';'T';'A';'Total'};
    qnm_450 = [qin_450;qr_450;qt_450;qa_450;qtot_450];
    qnm_500 = [qin_500;qr_500;qt_500;qa_500;qtot_500];
    qnm_550 = [qin_550;qr_550;qt_550;qa_550;qtot_550];
    qnm_800 = [qin_800;qr_800;qt_800;qa_800;qtot_800];
    qtotal = [qin_all;qr_all;qt_all;qa_all;qtot_all];
    Tablo2 = table(qnm_450,qnm_500,qnm_550,qnm_800,qtotal,'RowNames',qName)

    uitable('Data',Tablo2{:,:},'ColumnName',Tablo2.Properties.VariableNames,...
        'RowName',Tablo2.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

    % % saveas(fig11,'table.fig')
    % % saveas(fig11,'table.emf')
end
