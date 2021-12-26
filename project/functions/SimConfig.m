function SP = SimConfig(simP) 

%-------------------------------------------------------------------------%
%                               SimConfig.m
%-------------------------------------------------------------------------%
%   Function which sets the simulator parameters depending on the
%   configuration chosen from the user. A folder is created in the
%   current path where data and simulation's logfiles will be saved.
%   In the end, the consiguration is saved in the file SP.mat.
% 
%   INPUT: 
%    - simP           User input parameters (type: struct) 
%      .Mod_type      Modulation ( 1=QPSK , 2=16QAM , 3=64QAM )  
%      .CP_type       Cyclic prefix length ( 1=Normal , 2=Extended ) 
%      .T_sim         Simulation duration [s]
%      .Data_type     Flag for data bits generations ( 1=rand , 2=file )
%      .Sv_flag       Flag for data saving ( 1=ON , 0=OFF )
%      .TX_bw         Transmission Band ( 1.4, 3, 5, 10, 15, 20 MHz )
%      .P_shape       Impulse shape ( 1 = rect , 2 = raised cosine )   
%      .SNR           Signal to Noise Ratioe [dB]
%      .chan_type     Channel type ( 1=AWGN , 2=Multipath+AWGN )
%      .chan_ampl     Vetctor of channel's coefficients
%      .chan_tau      Vector of channel's latencies (tau)[s] 
%      .f_upsamp      Upsampling factor (simulation of continuos time)
% 
%   OUTPUT:
%    - SP             Simulation parameters (type: struct)
% 
%   Autore: Puccitelli Marco
%-------------------------------------------------------------------------%

    % Memorization of user input parameters
    % SP.Mod_type = simP.Mod_type;             
    % SP.CP_type = simP.CP_type;            	 
    % SP.T_sim = simP.T_sim;         	
    % SP.Data_type = simP.Data_type;  
    % SP.Sv_flag = simP.Sv_flag;
    % SP.TX_bw = simP.TX_bw;            
    % SP.P_shape = simP.P_shape;      
    % SP.f_upsamp = simP.f_upsamp;
    % SP.SNR = simP.SNR;
    % SP.chan_type = simP.chan_type;
    % SP.chan_ampl = simP.chan_ampl;
    % SP.chan_tau = simP.chan_tau;
    SP=simP;

    % Parametri del sistema
    SP.Ts = 1/(15e3*2048);             % Basic time unit [s]
    SP.delta_f = 15e3;                 % Spacing subcarriers [Hz]
    SP.N_sc_RB = 12;                   % Number of subcarriers for RB
    SP.T_frame = 0.5e-3;               % Subframe time duration[s]
    SP.N_frame = SP.T_sim/SP.T_frame;  % Total number of subframe per simulation
    SP.uplink = 1;

    switch SP.Mod_type
        case 1
            SP.modulazione = 'QPSK';
            SP.mod_bitsymb = 2;
        case 2
            SP.modulazione = '16QAM';
            SP.mod_bitsymb = 4;
        case 3
            SP.modulazione = '64QAM';
            SP.mod_bitsymb = 6;
        otherwise
            error('Valore Mod_type errato!');
    end    

    switch SP.TX_bw
        case 1.4e6
            SP.fs = 1.92e6;        % Sampling frequency
            SP.N_RB = 6;           % Number of occupied Resource Blocks
            SP.Msymb = 72;         % Number of occupied subcarriers
            SP.FFT_len = 128;      % Length of IFFT(tx)/FFT(rx)
        case 3e6
            SP.fs = 3.84e6;        % Sampling frequency
            SP.N_RB = 15;          % Number of occupied Resource Blocks
            SP.Msymb = 180;        % Number of occupied subcarriers
            SP.FFT_len = 256;      % Length of IFFT(tx)/FFT(rx)
        case 5e6
            SP.fs = 7.68e6;        % Sampling frequency
            SP.N_RB = 25;          % Number of occupied Resource Blocks
            SP.Msymb = 300;        % Number of occupied subcarriers
            SP.FFT_len = 512;      % Length of IFFT(tx)/FFT(rx)
        case 10e6
            SP.fs = 15.36e6;       % Sampling frequency
            SP.N_RB = 50;          % Number of occupied Resource Blocks
            SP.Msymb = 600;        % Number of occupied subcarriers
            SP.FFT_len = 1024;     % Length of IFFT(tx)/FFT(rx)
        case 15e6
            SP.fs = 23.04e6;       % Sampling frequency
            SP.N_RB = 75;          % Number of occupied Resource Blocks
            SP.Msymb = 900;        % Number of occupied subcarriers
            SP.FFT_len = 1536;     % Length of IFFT(tx)/FFT(rx)
        case 20e6
            SP.fs = 30.72e6;       % Sampling frequency
            SP.N_RB = 100;         % Number of occupied Resource Blocks
            SP.Msymb = 1200;       % Number of occupied subcarriers
            SP.FFT_len = 2048;     % Length of IFFT(tx)/FFT(rx)
        otherwise
            error('Valore di TX_bw errato!')
    end

    switch SP.CP_type
        case 1
            % Normal
            SP.CP = [160*SP.Ts 144*SP.Ts];      % Duration of cyclic prefix [s]
            SP.N_SCFDMAsymb = 14;               % SCFDMA symbols per subframe
            SP.Pilots=zeros(SP.Msymb,SP.N_SCFDMAsymb);
            SP.Flag_Pilots=zeros(SP.Msymb,SP.N_SCFDMAsymb);
            SP.Pilots(:,[4,11])=1+1i; 
            SP.Flag_Pilots(:,[4,11])=1; 
            SP.Reference_Symbols=[4,11];
        case 2
            % Extended
            SP.CP = 512*SP.Ts;                  % Duration of cyclic prefix [s]
            SP.N_SCFDMAsymb = 16;               % SCFDMA symbols per subframe
            SP.Pilots=zeros(SP.Msymb,SP.N_SCFDMAsymb);
            SP.Flag_Pilots=zeros(SP.Msymb,SP.N_SCFDMAsymb);
            SP.Flag_Pilots(:,[3,9])=1; 
            SP.Pilots(:,[3,9])=1+1i; 
            SP.Reference_Symbols=[3,9];
        otherwise
            error('Valore CP_type errato!');
    end
    % SP.Pilot_symbol=1+1i;

    % Folder creation for the simulator's outputs
    if isequal(SP.CHE,'LSi')
        SP.CHE_string_tag=[SP.CHE,num2str(SP.LSi_red_factor)];
    else
        SP.CHE_string_tag=[SP.CHE,num2str(SP.Version_CHE)];
    end
    SP.file_save = ['SCFDMA_',SP.modulazione,'_BW',num2str(SP.TX_bw),'_CP',num2str(SP.CP_type),'_',SP.SNR_comp,'_',SP.chan_type,num2str(SP.f_m),'_',SP.CHE_string_tag,'_NaRX',num2str(SP.n_ant_RX)];
    SP.OutputDir = ['.\results\'];       % directory file destination

    if (exist(SP.OutputDir,'dir') == 0)
        mkdir(SP.OutputDir);                   
    else
        file_save_corr=SP.file_save;
        while exist([SP.OutputDir,file_save_corr],'file')
            file_save_corr=[SP.file_save,'_s',num2str(floor(rand*100))];
        end
        SP.file_save=file_save_corr;
    end

    SP.noise_comp=1;

    % file_parametri = [SP.OutputDir,'SP_',SP.file_save,date,'.mat'];
    % save(file_parametri,'SP','-v6');

%End_Of_Function
end

%End_Of_File