%Program pemodelan inversi Data Transient Electromagnetic (TEM)/TDEM (Time-domain Electromagnetic) Central Loop Configuration
%menggunakan algoritma PSO (Particle Swarm Optimization)
%Mohammad Rheza Zamani
tic;
format long
clear all;
clc;
%Input data sintetik
t1=linspace(log10(10^-6),log(1),30);
t=10.^t1;
a=25;
I=1;
R = [100 1000 100];
thk = [100 100];
%Pehitungan data sintetik
TEM_sin =fwd_TEM(R,thk,t,a,I);
%Data Latihan
%Definisi ruang model 
npop = 100; %Jumlah dari model  
nlayer = 3; %Jumlah lapisan 
nitr = 1000; %Jumlah iterasi 
wmax = 0.9;
wmin = 0.5;
c1=2.05;
c2 =2.05;

%Batas atas dan bawah (Diatur 5 kali dari nilai parameter model)
%Batas bawah pencarian nilai resistivitas
LBR = [1 1 1];
%Batas atas pencarian nilai resistivitas
UBR = [2000 2000 2000];
%Batas bawah pencarian nilai ketebalan
LBT = [1 1];
%Batas atas pencarian nilai resistivitas
UBT = [200 200];
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

%Menhitung masing - masing dB/dt model semu dan misfit dari model
for ipop = 1 : npop 
    [TEM] = fwd_TEM(rho(ipop,:),thick(ipop,:),t,a,I);
    dBdt(ipop,:) = TEM;
    
    [misfit] = misfit_TEM(TEM_sin,dBdt(ipop,:));
    E(ipop) = misfit;
end
%Global best
idx = find(E ==min(E));
G_best_rho = rho(idx(1),:);
G_best_thick = thick(idx(1),:);
%Inversi
for itr = 1 : nitr
    w = wmin+((wmax-wmin)/nitr)*itr;
    for i = 1 : npop
        %Personal best / current best
        P_best_rho = rho;
        P_best_thick = thick;
        %Membuat komponen kecepatan
        %Rho
        for n1 = 1 : nlayer
            v_rho(1,n1) = w.*v_rho(i,n1) + c1.*rand.*(P_best_rho(n1) - rho(i,n1))+ c2.*rand.*(G_best_rho(n1) - rho(i,n1));
            rho_baru(1,n1) = rho(i,n1)+ v_rho(1,n1);
        if rho_baru(1,n1)<LBR(n1)
            rho_baru(1,n1) = LBR(n1);
        end
        if rho_baru(1,n1)>UBR(n1)
            rho_baru(1,n1) = UBR(n1);
        end
        end
        %Ketebalan
        for n2 = 1 : (nlayer-1)
            v_thk(1,n2) = w.*v_thk(i,n2) + c1.*rand.*(P_best_thick(n2) - thick(i,n2))+ c2.*rand.*(G_best_thick(n2) - thick(i,n2));
            thk_baru(1,n2) = thick(i,n2)+ v_thk(1,n2);
          if thk_baru(1,n2)<LBT(n2)
              thk_baru(1,n2) = LBT(n2);
          end
          if thk_baru(1,n2)>UBT(n2)
              thk_baru(1,n2) = UBT(n2);
          end
        end
        %Update pesonal best
        [dBdt_baru]= fwd_TEM(rho_baru,thk_baru,t,a,I);
        [err] = misfit_TEM(TEM_sin,dBdt_baru);
        if err<E(i)
            rho(i,:) = rho_baru(1,:);
            thick(i,:) = thk_baru(1,:);
            dBdt(i,:) = dBdt_baru(1,:);
            E(i) = err;
        end
    end
    Emin = 100;
     for ipop = 1 : npop
        if E(ipop)< Emin
            Emin = E(ipop);
            G_best_rho  = rho(ipop,:);
            G_best_thick  = thick(ipop,:);
            TEM_mod = dBdt(ipop,:);
        end
    end
    Egen(itr)=Emin;
end

    %Ploting model
r_plot = [0, R];
t_plot = [0,cumsum(thk),max(thk)*100];
r_mod = [0, G_best_rho];
Depth_mod = [0,cumsum(G_best_thick ),max(G_best_thick )*1000];

figure(1)
subplot(1,6,[1 3])
loglog(t,TEM_mod,'r',t,TEM_sin,'ob','MarkerSize',6,'LineWidth',2.5);
ylim([10^-20 10^-2])
legend({'Calculated Data','Observed Data'},'Color','none','FontWeight','Bold');
xlabel('Time (s)','FontSize',8,'FontWeight','Bold');
ylabel('dBdt (V/\itm^{2})','FontSize',8,'FontWeight','Bold');
title(['\bf \fontsize{10}\fontname{Times}TEM Respon  || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)]);

grid on
subplot(1,6,[5 6])
stairs(r_plot,t_plot,'--r','Linewidth',2);
hold on
stairs(r_mod,Depth_mod,'-b','Linewidth',2);
hold off
legend({'Synthetic Model','Inversion Model'},'Color','none','FontWeight','Bold','Location','southwest');
%legend({'Inversion Model'},'Color','none','FontWeight','Bold','Location','southwest');
axis([1 10^5 1 1000])
xlabel('Resistivity (\Omega.m)','FontSize',8,'FontWeight','Bold');
ylabel('Depth (m)','FontSize',8,'FontWeight','Bold');
title('\bf \fontsize{10} Model');
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
time = toc;
