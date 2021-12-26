%
% LTE_uplink_v01.m
% Simulator RUN Batch
%

%% Section for Transmissions, to change the settings see "Debug_interference.m"

N_sim = 1;     % Number of simulations
M_frame = 1;   % Number of frames

% Decomment if run only Run_LTE_uplink and viceversa
addpath('./functions')

% Simulate the real transmission, based on data of "Debug_interference.m"

BER_PER_GRAFICO={};
BER_PER_GRAFICO_CON_INTERFERENZA={};

% Initialization of struct BER_PER_GRAFICO_CON_INTERFERENZA
for i=1:3
    for j=1:numel(n_ant_RX_vect)
        n_ant_RX=n_ant_RX_vect(j);
        switch n_ant_RX
            case 1
            BER_PER_GRAFICO_CON_INTERFERENZA(i).ber1=0;  
            case 2
            BER_PER_GRAFICO_CON_INTERFERENZA(i).ber2=0;
            case 4
            BER_PER_GRAFICO_CON_INTERFERENZA(i).ber4=0;
            case 8
            BER_PER_GRAFICO_CON_INTERFERENZA(i).ber8=0;
        end
    end
end

% Combiner : chosen from input
combiner = input("Choose the combiner: 1 for MRC, 0 for MS\n");

if combiner
    combiner='MRC';
else
    combiner='MAX_SEL';
end

% Start of the simulation
for i=1:numel(n_ant_RX_vect)
    n_ant_RX=n_ant_RX_vect(i);
    for nue=1:numel(ue)
        if ue(nue).NNCellID
            fprintf('Utente %3.0f sta trasmettendo alla BS %5.0f \n',nue,ue(nue).NNCellID);
            EbTot=0;
            NoiseTot=0;
            indice_ebtot=1;
            
            % Choice of the interference: 1 exists, 0 not exist
            INTERFERENCE=1;
            LTE_uplink_v01;
            
            test=RES.BER;
            % Graph user per user of the BER
            switch n_ant_RX
                case 1
                    BER_PER_GRAFICO(nue).ber1_ant=mean(test,2);
                case 2
                    BER_PER_GRAFICO(nue).ber2_ant=mean(test,2);
                case 4
                    BER_PER_GRAFICO(nue).ber4_ant=mean(test,2);
                case 8
                    BER_PER_GRAFICO(nue).ber8_ant=mean(test,2);
            end
            
            % Struct of graphs depending on the number of interfering users
            for users_interfering=1:3
                if users_interfering==numel(ue(nue).u_int)+1
                    switch n_ant_RX
                        case 1
                        BER_PER_GRAFICO_CON_INTERFERENZA(users_interfering).ber1=BER_PER_GRAFICO_CON_INTERFERENZA(users_interfering).ber1...
                            +BER_PER_GRAFICO(nue).ber1_ant;      
                        case 2
                        BER_PER_GRAFICO_CON_INTERFERENZA(users_interfering).ber2=BER_PER_GRAFICO_CON_INTERFERENZA(users_interfering).ber2...
                            +BER_PER_GRAFICO(nue).ber2_ant;
                        case 4
                        BER_PER_GRAFICO_CON_INTERFERENZA(users_interfering).ber4=BER_PER_GRAFICO_CON_INTERFERENZA(users_interfering).ber4...
                            +BER_PER_GRAFICO(nue).ber4_ant;
                        case 8
                        BER_PER_GRAFICO_CON_INTERFERENZA(users_interfering).ber8=BER_PER_GRAFICO_CON_INTERFERENZA(users_interfering).ber8...
                            +BER_PER_GRAFICO(nue).ber8_ant;
                    end
                end
            end
        end
    end
end

for i=1:numel(indice)
    numero_utenti_interferenti=indice(i);
    for j=1:numel(n_ant_RX_vect)
        n_ant_RX=n_ant_RX_vect(j);
        switch n_ant_RX
            case 1
            BER_PER_GRAFICO_CON_INTERFERENZA(i).ber1=BER_PER_GRAFICO_CON_INTERFERENZA(i).ber1/numero_utenti_interferenti;  
            case 2
            BER_PER_GRAFICO_CON_INTERFERENZA(i).ber2=BER_PER_GRAFICO_CON_INTERFERENZA(i).ber2/numero_utenti_interferenti;
            case 4
            BER_PER_GRAFICO_CON_INTERFERENZA(i).ber4=BER_PER_GRAFICO_CON_INTERFERENZA(i).ber4/numero_utenti_interferenti;
            case 8
            BER_PER_GRAFICO_CON_INTERFERENZA(i).ber8=BER_PER_GRAFICO_CON_INTERFERENZA(i).ber8/numero_utenti_interferenti;
        end
    end
end

openvar ('BER_PER_GRAFICO_CON_INTERFERENZA');

