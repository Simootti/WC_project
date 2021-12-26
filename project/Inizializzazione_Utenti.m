%% Users initialization

% We must generate N UEs (User Equipment) which correspond to N IoT-devices
% associated with the 3 cells

% Nominal transmitting power Pr
Pt=10^(2.3-3); %23 dBm 
% Propagation coefficient Beta (variable from 2.5 to 6)
B=2.5;

% In general, you can generate N random numbers in the interval (a,b)
% with the formula r = a + (b-a).*rand(N,1)


% Inizialize the generator of random integers
for q=1:3*N
    rng('shuffle')
    rx(q)=x(2)+(x(3)-x(2)).*rand(1,1);
    ry(q)=y(1)+(y(2)-y(1)).*rand(1,1);
end

i=1;
%"j" index "Bs", "i" index "ue"
for j=1:3
    flag=1;
    while (flag)
        ue(i).user_id=i;
        ue(i).NBULSubcarrierSpacing = '15kHz'; % 3.75kHz, 15kHz
        ue(i).NNCellID = j;                    % Narrowband cell identity            
        ue(i).pos = [rx(i) ry(i)];
        ue(i).frameInterfering=[];
        ax=gca;
        figure(1)
        hold (ax,'on')
        scatter(ue(i).pos(1),ue(i).pos(2))
        i=i+1;
        if (i==N+1 || i==2*N+1 || i==3*N+1)
            flag=0;
        end
    end
end

for j=1:3
    Bs(j).Number_of_users=0;
    for i=1:numel(ue)
        ue(i).Pt=Pt;
        Bs(j).Number_of_users=0; 
        % X represents the distances between the user ue(i)
        % which belongs to the first cell and to the BS
        X = [Bs(j).pos;ue(i).pos];
        Bs(j).d(i) = pdist(X,'euclidean');
        Bs(j).Pr_avg(i)=0;
        Bs(j).Pn_avg(i)=0;
        % For each "ue" we memorize the distance, the pathloss and the
        % received power
        Bs(j).PL(i)= PathLoss(Bs(j).d(i),f_c,B);
        % In the Friis formula, we need the G of ant Tx and Rx
        % We considered the antennas' gains Gtx=0dB (1) and
        % Grx=17 dB (50.12 dB) and NF(Noise Figure)=3dB
        Gtx=10^(0/10);
        Grx=10^(1.7);
        NF=10^.3;
        
        switch shadowing
            
            case 'uniforme'
                ue(i).sh=10^(shadowing_uniforme(std_db)/10);
                
            case 'non_uniforme'
                ue(i).sh=10^(shadowing_non_uniforme(ue(i).pos(1),ue(i).pos(2),std_db,shadow)/10);
        end

        Bs(j).Pr_nominale(i) = 10*log10((ue(i).Pt*Gtx*Grx*1e3)/(Bs(j).PL(i)*NF*ue(i).sh));
        Bs(j).Pr_avg(i)=0;
        Bs(j).Pn_avg(i)=0;          
    end   
end