%Program pemodelan inversi Data Transient Electromagnetic (TEM)/TDEM (Time-domain Electromagnetic) Central Loop Configuration
%menggunakan algoritma FPA (Flower Pollination Algorithm)
%Mohammad Rheza Zamani
tic
clear all;
clc;
%Input data sintetik
t1=linspace(log10(10^-6),log(1),30);
t=10.^t1;
a=25;
I=1;
R = [100 1000 10];
thk = [100 100];
%Pehitungan data sintetik
TEM_sin =fwd_TEM(R,thk,t,a,I);
%Data Latihan
%t = [ 0.00004007 0.0001011 0.0001621 0.0002378 0.0003601 0.000527 0.0008145 0.001558];
%TEM_sin = [ 0.0000054811749 0.00000098437 0.000000396178 0.00000018854 0.000000083357 0.000000038893 0.0000000156912 0.0000000034128];
%Definisi ruang model 
npop = 100; %Jumlah dari model  
nlayer = 3; %Jumlah lapisan 
nitr = 200; %Jumlah iterasi 
prp=0.8;                                %probability
%Batas atas dan bawah (Diatur 5 kali dari nilai parameter model)
%Batas bawah pencarian nilai resistivitas
LBR = [1 1 1];
%Batas atas pencarian nilai resistivitas
UBR = [300 3000 30];
%Batas bawah pencarian nilai ketebalan
LBT = [1 1];
%Batas atas pencarian nilai resistivitas
UBT = [300 300];
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
idx = find(E ==min(E));
rho_best = rho(idx(1),:);
thick_best = thick(idx(1),:);
%Inversi
for iitr =  1 : nitr
    for i = 1 : npop
        if rand < prp
            n = 1;
            m1 = nlayer;
            beta = 1.5;
            [Lr]=levy(n,m1,beta);%levy distribution
            for n = 1 : nlayer
                npopr(1,n)=rho(i,n)+(0.1.*Lr(n).*(rho(i,n)-rho_best(n)));
                if npopr(1,n) < LBR(n);
                     npopr(1,n) = LBR(n);
                end
                if npopr(1,n) > UBR(n);
                     npopr(1,n) = UBR(n);
                end
            end
            m2 = nlayer -1;
            [Lt]=levy(n,m2,beta);%levy distribution
            for n = 1 : (nlayer-1)
                npopt(1,n)=thick(i,n)+(0.1.*Lt(n).*(thick(i,n)-thick_best(n)));
                if npopt(1,n) < LBT(n);
                     npopt(1,n) = LBT(n);
                end
                if npopt(1,n) > UBT(n);
                    npopt(1,n) = UBT(n);
                 end
            end
          else
            epsilon=rand;
            JK=randperm(npop);
            for n = 1 : nlayer
                npopr(1,n)=rho(i,n)+epsilon*(rho(JK(1),n)-rho(JK(2),n));
                 if npopr(1,n) < LBR(n);
                     npopr(1,n) = LBR(n);
                 end
                 if npopr(1,n) > UBR(n);
                    npopr(1,n) = UBR(n);
                 end
            end
            for n = 1 : (nlayer-1)
                npopt(1,n)=thick(i,n)+epsilon*(thick(JK(1),n)-thick(JK(2),n));
                 if npopt(1,n) < LBT(n);
                     npopt(1,n) = LBT(n);
                 end
                 if npopt(1,n) > UBT(n);
                     npopt(1,n) = UBT(n);
                end
            end
        end
        rho_new=npopr;
        thick_new=npopt;
        [dBdt_new]= fwd_TEM(rho_new,thick_new,t,a,I);
        [err] = misfit_TEM(TEM_sin,dBdt_new);
        %Update model dan error jika lebih baik
        if err<E(i)
            rho(i,:) = rho_new(1,:);
            thick(i,:) = thick_new(1,:);
            dBdt(i,:) = dBdt_new(1,:);
            E(i) = err;
        end
    end
    Emin = 100;
    for ipop = 1 : npop
        if E(ipop)< Emin
            Emin = E(ipop);
            rho_best = rho(ipop,:);
            thick_best = thick(ipop,:);
            TEM_mod = dBdt(ipop,:);
        end
    end
    Egen(iitr)=Emin;
end

    %Persiapan ploting
%Ploting model
r_plot = [0, R];
t_plot = [0,cumsum(thk),max(thk)*100];
r_mod = [0,rho_best];
Depth_mod = [0,cumsum(thick_best),max(thick_best)*100];

figure(1)
subplot(1,6,[1 3])
loglog(t,TEM_mod,'r',t,TEM_sin,'ob','MarkerSize',6,'LineWidth',2.5);
ylim([10^-20 10^-2])
legend({'Calculated Data','Observed Data'},'Color','none','FontWeight','Bold');
xlabel('AB/2 (m)','FontSize',8,'FontWeight','Bold');
ylabel('dBdt (V/\itm^{2}) (Ohm.m)','FontSize',8,'FontWeight','Bold');
title(['\bf \fontsize{10}\fontname{Times}TEM Respon  || Misfit : ', num2str(Egen(iitr)),' || iteration : ', num2str(iitr)]);
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
title('\bf \fontsize{12} Misfit Graph ');
grid on
time = toc;

%Levy Function
function [z] = levy(n,m,beta)
    num = gamma(1+beta)*sin(pi*beta/2);
    
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2);

    sigma_u = (num/den)^(1/beta);

    u = normrnd(0,sigma_u^2,n,m); 
    
    v = normrnd(0,1,n,m);

    z = u./(abs(v).^(1/beta));
end