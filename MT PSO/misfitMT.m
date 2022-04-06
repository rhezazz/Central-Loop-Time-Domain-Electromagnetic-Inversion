function [misfit]=misfitMT(rho_obs,phase_obs,rho_cal,phase_cal)
l=length(rho_obs);
d2r=pi/180;
lamda = 0.1;
for i=1:l
    m(i)=(abs(log10(rho_cal(i)/rho_obs(i)))+abs(d2r*phase_cal(i)-d2r*phase_obs(i)));
     %m(i)=(abs(log10(rho_cal(i)/rho_obs(i)))+lamda*(abs(d2r*phase_cal(i)-d2r*phase_obs(i))));
end
misfit=sum(m)/(l);
end