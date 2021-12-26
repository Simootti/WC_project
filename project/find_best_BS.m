%% Function which associates users to a Base Station

max=Sensitivity;
flag=0;
for n_bs=1:3 % Number of Base Stations (BS)
    if (Bs(n_bs).Pr_nominale(nue)>=max && Bs(n_bs).d(nue)<=R)
        max=Bs(n_bs).Pr_nominale(nue);
        flag=1;
        nbs_pref=n_bs;
        Bs(n_bs).Number_of_users=Bs(n_bs).Number_of_users + 1;
    end
end

if flag
    ue(nue).NNCellID=nbs_pref;
else
    ue(nue).NNCellID=0;
end

clear max nbs_pref flag
return