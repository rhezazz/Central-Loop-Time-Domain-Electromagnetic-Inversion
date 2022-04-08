    %Program pemodelan inversi Data Transient Electromagnetic (TEM)/TDEM (Time-domain Electromagnetic) Central Loop Configuration
%menggunakan algoritma Fireflies Algorithm (FA)
%Mohammad Rheza Zamani
tic;
format bank
clear all;
clc;
%Input data sintetik
t1=linspace(log10(10^-6),log(1),30);
t=10.^t1;
a=25;
I=1;
R = [100 500 100];
thk = [100 500];
%Pehitungan data sintetik
TEM_sin =fwd_TEM(R,thk,t,a,I);
%Data Latihan
%t = [ 0.00004007 0.0001011 0.0001621 0.0002378 0.0003601 0.000527 0.0008145 0.001558];
%TEM_sin = [ 0.0000054811749 0.00000098437 0.000000396178 0.00000018854 0.000000083357 0.000000038893 0.0000000156912 0.0000000034128];
%Definisi ruang model 
npop = 100; %Jumlah dari model  
nlayer = 3; %Jumlah lapisan 
nitr = 500; %Jumlah iterasi 
%Batas atas dan bawah (Diatur 3 kali dari nilai parameter model)
%Batas bawah pencarian nilai resistivitas
LBR = [1 1 1];
%Batas atas pencarian nilai resistivitas
UBR = [500 2500 500];
%Batas bawah pencarian nilai ketebalan
LBT = [1 1];
%Batas atas pencarian nilai resistivitas
UBT = [500 2500];
alpha = 0.2;
betha0 = 1;
gamma = 0.8;
damp = 0.99;
%Membuat model awal acak
for ipop = 1 : npop
    rho(ipop , :) = LBR + rand*(UBR - LBR);
    thick(ipop, :) = LBT + rand*(UBT - LBT);
end
%Menhitung masing - masing dB/dt model semu dan misfit dari model
for ipop = 1 : npop 
    [TEM] = fwd_TEM(rho(ipop,:),thick(ipop,:),t,a,I);
    dBdt(ipop,:) = TEM;
    
    [misfit] = misfit_TEM(TEM_sin,dBdt(ipop,:));
    E(ipop) = misfit;
end

%Inversi
for itr = 1 : nitr
    for i =  1 : npop
        j = randi(npop,1);
        while i == j
            j = randi(npop,1);
        end
        if E(i)<E(j)
            dr = norm((rho(i,:)-rho(j,:)));
            dt = norm((thick(i,:)-thick(j,:)));
            %Calculated new model with determine a new position
            %Random vector position for resistivity model
            for n1 = 1 : nlayer
                rho_baru(1,n1) = rho(i,n1) +betha0.*exp(-gamma*(dr)^2).*(rho(j,n1)-rho(i,n1))+ (alpha*(rand-0.5)*abs((UBR(n1)-LBR(n1))));
                if rho_baru(1,n1) < LBR(n1);
                     rho_baru(1,n1) = LBR(n1);
                end
                if rho_baru(1,n1) > UBR(n1);
                     rho_baru(1,n1) = UBR(n1);
                end
            end
            %Random vector position for thick model
            
            for n2 = 1 : (nlayer-1)
                thk_baru(1,n2) = thick(i,n2) +betha0.*exp(-gamma*(dt).^2).*(thick(j,n2)-thick(i,n2))+ (alpha*(rand-0.5)*abs((UBT(n2)-LBT(n2))));
                if thk_baru(1,n2) < LBT(n2);
                     thk_baru(1,n2) = LBT(n2);
                end
                if thk_baru(1,n2) > UBT(n2);
                    thk_baru(1,n2) = UBT(n2);
                end
            end 
        else
            rho_baru(1,:) = rho(i,:);
            thk_baru(1,:) = thick(i,:);
        end
        [dBdt_baru]= fwd_TEM(rho_baru,thk_baru,t,a,I);
        [err] = misfit_TEM(TEM_sin,dBdt_baru);
        %Update model dan error jika lebih baik
        if err<E(i)
            rho(i,:) = rho_baru(1,:);
            thick(i,:) = thk_baru(1,:);
            dBdt(i,:) = dBdt_baru(1,:);
            E(i) = err;
        end
      
    end
     Emin = 1000;
    for ipop = 1 : npop
        if E(ipop)< Emin
            Emin = E(ipop);
            rho_model = rho(ipop,:);
            thk_model = thick(ipop,:);
            TEM_mod = dBdt(ipop,:);
        end
    end
    Egen(itr)=Emin;
    alpha = alpha*damp;
end
time = toc
        %Ploting model
r_plot = [0, R];
t_plot = [0,cumsum(thk),max(thk)*100];
r_mod = [0,rho_model];
Depth_mod = [0,cumsum(thk_model),max(thk_model)*100];

figure(1)
subplot(1,6,[1 3])
loglog(t,TEM_mod,'r',t,TEM_sin,'ob','MarkerSize',6,'LineWidth',2.5);
ylim([10^-20 10^-2])
legend({'Calculated Data','Observed Data'},'Color','none','FontWeight','Bold');
xlabel('Time (ms)','FontSize',8,'FontWeight','Bold');
ylabel('dBdt (V/\itm^{2}) (Ohm.m)','FontSize',8,'FontWeight','Bold');
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
xlabel('Resistivity (Ohm.m)','FontSize',8,'FontWeight','Bold');
ylabel('Depth (m)','FontSize',8,'FontWeight','Bold');
title('\bf \fontsize{10} Model');
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
