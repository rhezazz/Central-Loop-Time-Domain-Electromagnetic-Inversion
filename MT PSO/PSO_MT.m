%Program pemodelan inversi kurva sounding MT 1-D dengan
%menggunakan algoritma Particle Swarm Optimization
%Mohammad Rheza Zamani
tic
clear all;
clc;
%Data sintetik
R = [500 100 1000];
thk = [500 1500];
freq = logspace(-3,3,50);
T = 1./freq;
[app_sin, phase_sin] = modelMT(R, thk ,T);

%Definisi ruang model 
npop = 100; %Jumlah dari model  
nlayer = 3; %Jumlah lapisan 
nitr = 1000; %Jumlah iterasi 
%Ruang pencarian diatur sebesar 5 kali dari model data  sintetik
%Batas bawah pencarian nilai resistivitas
LBR = [1 1 1];
%Batas atas pencarian nilai resistivitas
UBR = [2000 2000 2000];
%Batas bawah pencarian nilai ketebalan
LBT = [1 1];
%Batas atas pencarian nilai resistivitas
UBT = [2000 2000];
%Paramter inversi
wmax = 0.9;
wmin = 0.5;
c1=2.05;
c2 =2.05;
%Membuat model awal acak
for ipop = 1 : npop
    rho(ipop , :) = LBR + rand*(UBR - LBR);
    thick(ipop, :) = LBT + rand*(UBT - LBT);
end
%Hitung velocity awal dan posisi awal
for ipop = 1 : npop
    for imod =  1 : nlayer
        %v_rho(ipop,imod) = 0.5.*(min(rho(ipop,:)) + rand*(max(rho(ipop,:)) - min(rho(ipop,:))));
        v_rho(ipop,imod) = 0;
    end
    for imod = 1 : nlayer -1
        %v_thk(ipop,imod) = 0.5.*(min(thick(ipop,:))) + rand*(max(thick(ipop,:)) - min(thick(ipop,:)));
        v_thk(ipop,imod) = 0;
    end
end
for ipop=1:npop
    [apparentResistivity, phase_baru]=modelMT(rho(ipop,:),thick(ipop,:),T);
     app_mod(ipop,:)=apparentResistivity;
     phase_mod(ipop,:)=phase_baru;
     
    [misfit]=misfitMT(app_sin,phase_sin,app_mod(ipop,:),phase_mod(ipop,:));
    E(ipop)=misfit;
end
%Global best
idx = find(E ==min(E));
G_best_rho = rho(idx(1),:);
G_best_thick = thick(idx(1),:);
%Inversi
for itr = 1 : nitr
    w = wmax-((wmax-wmin)/nitr)*itr;
    for i = 1 : npop
        P_best_rho = rho;
        P_best_thick = thick;
        %Membuat komponen kecepatan
        %Rho
        for imod = 1 : nlayer
            v_rho(1,imod) = w.*v_rho(i,imod) + c1.*rand.*(P_best_rho(i,imod) - rho(i,imod))+ c2.*rand.*(G_best_rho(imod) - rho(i,imod));
            rho_baru(1,imod) = rho(i,imod)+ v_rho(1,imod);
        if rho_baru(1,imod)<LBR(imod)
            rho_baru(1,imod) = LBR(imod);
        end
        if rho_baru(1,imod)>UBR(imod)
            rho_baru(1,imod) = UBR(imod);
        end
        end
        %Ketebalan
        for imod = 1 : (nlayer-1)
            v_thk(1,imod) = w.*v_thk(i,imod) + c1.*rand.*(P_best_thick(i,imod) - thick(i,imod))+ c2.*rand.*(G_best_thick(imod) - thick(i,imod));
            thick_baru(1,imod) = thick(i,imod)+ v_thk(1,imod);
          if thick_baru(1,imod)<LBT(imod)
              thick_baru(1,imod) = LBT(imod);
          end
          if thick_baru(1,imod)>UBT(imod)
              thick_baru(1,imod) = UBT(imod);
          end
        end
        %Update pesonal best
        [apparentResistivity_baru, phase_baru]=modelMT(rho_baru,thick_baru,T);
        app_mod_baru = apparentResistivity_baru;
        phase_mod_baru(i,:) = phase_baru;
        [E_baru] = misfitMT(app_sin,phase_sin,app_mod_baru, phase_mod_baru);
        if E_baru<E(i)
            rho(i,:) = rho_baru(1,:);
            thick(i,:) = thick_baru(1,:);
            app_mod(i,:) = app_mod_baru;
            phase_mod(i,:) = phase_mod_baru(1,:);
            E(i) = E_baru;
        end
    end
    Emin = 100;
     for ipop = 1 : npop
        if E(ipop)< Emin
            Emin = E(ipop);
            G_best_rho  = P_best_rho(ipop,:);
            G_best_thick  = P_best_thick(ipop,:);
            app_model= app_mod(ipop,:);
            phase_model = phase_mod(ipop,:);
        end
    end
    Egen(itr)=Emin;
end
toc
    %Persiapan Ploting
    rho_plot = [0 R];
    thk_plot = [0 cumsum(thk) max(thk)*10000];
    rhomod_plot = [0 G_best_rho];
    thkmod_plot = [0 cumsum(G_best_thick ) max(G_best_thick )*10000];
    %Plotting
    figure(1)
    subplot(2, 2, 1)
    loglog(T,app_sin,'.b',T,app_model,'r','MarkerSize',12,'LineWidth',1.5);
    axis([10^-3 10^3 1 10^4]);
    legend({'Synthetic Data','Calculated Data'},'EdgeColor','none','Color','none','FontWeight','Bold','location','southeast');
    xlabel('Periods (s)','FontSize',12,'FontWeight','Bold');
    ylabel('App. Resistivity (Ohm.m)','FontSize',12,'FontWeight','Bold');
    title(['\bf \fontsize{10}\fontname{Times}Period (s) vs Apparent Resistivity (ohm.m)  || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)]);
    grid on
    
    subplot(2, 2, 3)
    loglog(T,phase_sin,'.b',T,phase_model,'r','MarkerSize',12,'LineWidth',1.5);
    axis([10^-3 10^3 0 90]);
    set(gca, 'YScale', 'linear');
    legend({'Synthetic Data','Calculated Data'},'EdgeColor','none','Color','none','FontWeight','Bold');
    xlabel('Periods (s)','FontSize',12,'FontWeight','Bold');
    ylabel('Phase (deg)','FontSize',12,'FontWeight','Bold');
    title(['\bf \fontsize{10}\fontname{Times}Period (s) vs Phase (deg)  || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)]);
    grid on
    
    subplot(2, 2, [2 4])
    stairs(rho_plot,thk_plot,'--b','Linewidth',1.5);
    hold on
    stairs(rhomod_plot ,thkmod_plot,'-r','Linewidth',2);
    hold off
     legend({'Synthetic Model','Calculated Model'},'EdgeColor','none','Color','none','FontWeight','Bold','Location','SouthWest');
    axis([1 10^4 0 5000]);
    xlabel('Resistivity (Ohm.m)','FontSize',12,'FontWeight','Bold');
    ylabel('Depth (m)','FontSize',12,'FontWeight','Bold');
    title(['\bf \fontsize{10}\fontname{Times}Model']);
    subtitle(['\rho_{1} = ',num2str(G_best_rho(1)),' || \rho_{2} = ',num2str(G_best_rho(2)),' || \rho_{3} = ',num2str(G_best_rho(3)),' || thick_{1} = ',num2str(G_best_thick(1)),' || thick_{2} = ',num2str(G_best_thick(2))],'FontWeight','bold')
    set(gca,'YDir','Reverse');
    set(gca, 'XScale', 'log');
    set(gcf, 'Position', get(0, 'Screensize'));
    grid on

%plot misfit
figure(2)
plot(1:nitr,Egen,'r','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('RSME','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Misfit ');
grid on