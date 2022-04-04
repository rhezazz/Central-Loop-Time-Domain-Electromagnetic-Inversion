%Fungsi Objektif Pemodelan Inversi Data TEM/TDEM
function [misfit] = misfit_TEM(TEM_obs, TEM_cal)
    ls = length(TEM_obs);
    for j = 1 : ls
        m(j) = ((TEM_cal(j) - TEM_obs(j))/TEM_obs(j))^2;
    end
    misfit = sqrt((1/ls)*sum(m));
end
