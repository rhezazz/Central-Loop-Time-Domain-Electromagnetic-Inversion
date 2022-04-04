%Forward Modelling TEM Data (Transient Electromagnetic Central Loop)
%Mohammad Rheza Zamani
%Reference : Christensen, N. B. A Generic 1-D Imaging Method for Transient Electromagnetic Data. Geophysics. 2002. 67, 438-447.
function [dBdt] = fwd_TEM(rho,thk,t,a,I)
u0=4*pi*10^(-7); 
alpha = 0.6;
c = 1.2;
con = 1./rho;
z = [0 cumsum(thk) inf];
for i = 1 : length(t)
    av_con(i) = mean(con);
end

for k = 1 : 10
    d = sqrt((c*t)./(av_con*u0));
    F1 = zeros(length(t),(length(z)-1));
    F2 = zeros(length(t),(length(z)-1));
    F=zeros(length(t),length(z)-1);
    for i = 1 : length(t)
        for j = 1 : (length(z)-1)
            if z(j) <= d(i) 
                F1(i,j) = (2-(z(j)/d(i))).*(z(j)/d(i));
            elseif z(j)> d(i)
                F1(i,j) = 1;
            end
           if z(j+1) <= d(i) 
                F2(i,j) = (2-(z(j+1)/d(i))).*(z(j+1)/d(i));
            elseif z(j+1)> d(i)
                F2(i,j) = 1;
            end
            F(i,j)=F2(i,j)-F1(i,j);    
        end
    end
                app_con1 = con'.*F';
                app_con = sum(app_con1);
                av_con = alpha*app_con + (1-alpha)*av_con;
                
end
tetha = sqrt((u0*av_con)./(4*t));
%Faraday's Law dB/dt = -u0*(dH/dt)
p = (I./(av_con*(a^3)));
q = 3*erf(tetha*a);
r = (2/sqrt(pi)).*tetha.*a.*(3+2.*((tetha*a).^2));
s = exp(-((tetha*a).^2));
dBdt = p.*(q - r.*s);
end