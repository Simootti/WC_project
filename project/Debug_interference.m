%% Debug Interference

Mod_type = 1;                % 1=QPSK , 2=16QAM , 3=64QAM   
N_sim = 1;                   % Number of simulations
To=290;                      % Temperature [°Kelvin]
K = physconst('Boltzmann');  % Boltzmann's constant
TX_bw = 20e6;                % 1.4 , 3 , 5 , 10 , 15 , 20 [MHz] 
SNR = linspace(0,15,6);      % Signal to Noise Ratio [dB] for the graph
SNR_comp = '3GPP' ;          % SNR computation: 'MEAS', '3GPP'
chan_type ='RAYL';           % 'AWGN','RAYL','EPA','EVA','ETU'

%Channel Estimation
CHE = 'LSf';                 % 'IDEAL','LS','MMSE'
LSi_red_factor=4;
NVLE = 'IDEAL';
Version_CHE=1;               % 0, 1 (ideal L), 
f_m = 5;                     % Max Doppler (5 Hz)    

% Before doing the "real" simulation, it must be evaluated the BS to which
% the user will connect, i.e. the one with 
% the higher received power from the "ue"

for nue=1:numel(ue)
    find_best_BS;
end

% Finds all the interferences that each user could have during  
% the transmission
for nue=1:numel(ue)
    ue(nue).number_users_interfering=0;
    ue(nue).users_interfering=[];
    ue(nue).flag_tx_int=0;
    k=1;
    if(ue(nue).NNCellID)
        fprintf('Utente %5.0f dovrà trasmettere alla BS %5.0f \n',nue,ue(nue).NNCellID);
        for user_int=1:numel(ue)
            if(ue(user_int).NNCellID)
                X = [ue(nue).pos;ue(user_int).pos];
                if(ue(user_int).NNCellID ~= ue(nue).NNCellID)
                    ue(nue).dist_wrto_user(user_int) = pdist(X,'euclidean');
                    fprintf("Utente %3.0f della Bs %3.0f interferirà con l'utente %3.0f\n",...
                        user_int,ue(user_int).NNCellID,nue);
                    ue(nue).number_users_interfering=ue(nue).number_users_interfering + 1;
                    ue(nue).users_interfering(k)=user_int;
                    k=k+1;
                end
            end
        end
    else
        fprintf('Utente %5.0f non può trasmettere \n',nue);
    end
        
end

indice=zeros(1,3);

% Random choice of interfering users
for nue=1:numel(ue)
    if ue(nue).number_users_interfering
        if round(rand())
            primo_utente_casuale=randperm(length(ue(nue).users_interfering),1);
            ue(nue).u_int(1)=ue(nue).users_interfering(primo_utente_casuale);
            if round(rand())
                indice_secondo_utente_casuale=randperm(length(ue(nue).users_interfering),1);
                secondo_utente_casuale=ue(nue).users_interfering(indice_secondo_utente_casuale);
                while ue(secondo_utente_casuale).NNCellID == ue(ue(nue).u_int(1)).NNCellID
                    indice_secondo_utente_casuale=randperm(length(ue(nue).users_interfering),1);
                    secondo_utente_casuale=ue(nue).users_interfering(indice_secondo_utente_casuale);
                end
                ue(nue).u_int(2)=secondo_utente_casuale;
                ue(nue).number_users_interfering=2;
                indice(3)=indice(3)+1;
            else
                ue(nue).number_users_interfering=1;
                indice(2)=indice(2)+1;
            end
        else
            ue(nue).number_users_interfering=0;
            indice(1)=indice(1)+1;
        end
    end
end

fprintf("%s utenti con 0 utenti interferenti\n%s utenti con 1 utente interferente\n%s utenti con 2 utenti interferenti\n",...
    num2str(indice(1)),num2str(indice(2)),num2str(indice(3)));


n_ant_RX=8;

wtb = waitbar(0,'Attendere la creazione dei frame interferenti...');

% For loop which orders the users' transmissions from the one with more
% interference the one with less interference
for nue=1:numel(ue)
     for index=1:ue(nue).number_users_interfering  %For each user
        nue_int=ue(nue).u_int(index);  %In the vector of interfering users, i create an interf. frame
        if ue(nue_int).flag_tx_int %If that user has already tx interf. frame, skip interf. user
        else
            Generating_interfering_Frame;
            waitbar((nue)/(numel(ue)),wtb,['Utente interferente (', num2str(index),...
                    '/',num2str(ue(nue).number_users_interfering),'),',...
                        'Utente (',num2str(nue),'/',num2str(numel(ue)),') ',num2str(ceil((nue)/(numel(ue))*100)) '%']);
        end
     end
end

close(wtb);

% In order to see the most interfered users, could be useful a struct ordering 
% wrto the number of interfering users
ue_int_order=nestedSortStruct(ue,'number_users_interfering',-1);
openvar('ue')
