%% Generating Interfering Frame

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INPUT: 
%    - simP           User input parameters (type: struct) 
%      .Mod_type      Modulation ( 1=QPSK , 2=16QAM , 3=64QAM )  
%      .CP_type       Cyclic prefix length ( 1=Normal , 2=Extended ) 
%      .T_sim         Simulation duration [s]
%      .Data_type     Flag for data bits generation ( 1=rand , 2=file )
%      .Sv_flag       Flag for saving data ( 1=ON , 0=OFF )
%      .TX_bw         Transmission Band ( 1.4, 3, 5, 10, 15, 20 MHz )
%      .P_shape       Impulse shape ( 1=rect , 2=coseno rialzato )   
%      .SNR           Signal to Nois Ratio [dB]
%      .chan_type     Channel Type
%      .chan_ampl     Vector of channel coefficients
%      .chan_tau      Vector of channel latencies (tau) [s] 
%      .f_upsamp      Upsampling factor (simulates the continuos time)
%
%      .n_ant_RX      Number of BS antennas
%                   
% 
%   OUTPUT:
%    - PAPR modulation
%    - BER for modulation and SNR selected
%    - SER for modulation and SNR selected
%    - SNR measured at the receiver 
%    - Costellation (I,Q)
% 
%   NB: the script uses some of the following functions  
%       (that must be in the same folder of the simulator):
%       
%       - SimConfig.m       Returns the configuration parameters 
%                           of the simulation
%       - DataGen.m         Returns the bits that must be transmitted
%       - ConstMap.m        Mapping of bit/symbols for the chosen
%                           modulation
%       - AddCP.m           Adds the cyclic prefix      
%       - RemCP.m           Removes the cyclic prefix       
%       - Decisor.m         Estimates the received symbol
%       - ConstDemap.m      Demapping of received symbols in bit
%       - SimLog.m          Generates the log file of the simulation
%       - constIQplot.m     Represents the transmitted symbols
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc, clear all, close all

% from Run_...
% Mod_type
% N_sim
% TX_bw 
% SNRV 
% SNR_comp
% chan_type 
% CHE
% f_m
% n_ant_RX

%-------------------------------------------------------------------------%
%                          Simulator Parameters 
%-------------------------------------------------------------------------%

% Parameters Modulator/Demodulator

if not(exist('Mod_type'))
    simP.Mod_type = 1;               % 1=QPSK , 2=16QAM , 3=64QAM   
else
    simP.Mod_type = Mod_type;        % 1=QPSK , 2=16QAM , 3=64QAM   
end
simP.CP_type = 1;                    % Cyclic Prefix 1=Normal , 2=Extended
simP.T_sim = 0.01;                   % Simulation Duration [s]
if not(exist('N_sim'))               % Number of simulations
    simP.N_sim = 10;
else
    simP.N_sim = N_sim;
end
simP.Data_type = 1;                  % Flag of data bit ( 1=rand , 2=file )
if not(exist('TX_bw'))
    simP.TX_bw = 10;                 % 1.4 , 3 , 5 , 10 , 15 , 20 [MHz] 
else
    simP.TX_bw = TX_bw;
end
simP.P_shape = 1;                    % 1=rect , 2=raised cosine
simP.f_upsamp = 4;                   % Upsampling factor

% Parametri Canale
if not(exist('SNR'))
    simP.SNR = [30];                 % SNR [dB]
else
    simP.SNR = SNR;
end
if not(exist('SNR_comp'))
    simP.SNR_comp = '3GPP';          % SNR computation: 'MEAS', '3GPP'
else
    simP.SNR_comp = SNR_comp;
end
if not(exist('chan_type'))           % Rayleigh
    simP.chan_type = 'RAYL';         % 'AWGN','RAYL','EPA','EVA','ETU'
else
    simP.chan_type = chan_type;
end
if not(exist('CHE'))                 % Channel Estimation->Least Square
    simP.CHE = 'LS';                 % 'IDEAL','LS','LSf','MMSE', 'LSi'
else
    simP.CHE = CHE;
end
if not(exist('NVLE'))                % (Noise)
    simP.NVLE = 'IDEAL';             % 'IDEAL','LS','LSf','MMSE'
else
    simP.NVLE = NVLE;
end
if not(exist('Version_CHE'))
    simP.Version_CHE = 0;            % 0, 1 (ideal L)
else
    simP.Version_CHE = Version_CHE;
end
if not(strcmp(simP.CHE,'IDEAL'))&&not(strcmp(simP.NVLE,'IDEAL'))
    simP.noise_comp=1;
else
    simP.noise_comp=0;
end
if not(exist('LSi_red_factor'))
    simP.LSi_red_factor = 20;         % 0, 1 (ideal L)
else
    simP.LSi_red_factor = LSi_red_factor;
end

simP.chan_ampl = [1 .2];              % Channel coefficients vector (2 rays)
simP.chan_tau = [0 5e-6];             % Vector of channel latencies (tau)[s]
simP.f0 = 2e+9;                       % carrier frequency, 2GHz
if not(exist('f_m'))
    simP.f_m = 5;                     % Max Doppler
else
    simP.f_m = f_m;
end

%% Parametri BS
if not(exist('n_ant_RX'))
     simP.n_ant_RX = 1;
else
    simP.n_ant_RX = n_ant_RX;
end

% Parametri simulazione
simP.Sv_flag = 1;       % Flag data saving ( 1=ON , 0=OFF, 2=ON+DATA )
simP.plot_fig = 0;

%-------------------------------------------------------------------------%
%                   Simulator configuration SCFDMA
%-------------------------------------------------------------------------%
SP = SimConfig(simP); 

% Inizialization of counters e variables input/output
clear RES
RES.SEcount = zeros(length(SP.SNR),SP.N_sim);
RES.BEcount = zeros(length(SP.SNR),SP.N_sim);
RES.Scount = zeros(length(SP.SNR),SP.N_sim);
RES.Bcount = zeros(length(SP.SNR),SP.N_sim);
RES.PAPR = zeros(length(SP.SNR),SP.N_sim);
RES.Err_che = zeros(length(SP.SNR),SP.Msymb);
Input.DATA = [];
Input.SY_TX = [];
Output.SY_RX = [];
Output.SY_stima = [];
Output.SY_err = [];
Output.DATA = [];
snr_misurato = zeros(1,N_sim);
snr_misurato_symb = zeros(1,SP.N_SCFDMAsymb);
SGNL={};
SGNL_nn={};
% SNR is calculated for each symbol and for each frame, so all the elements
% are summed in the end and the SNR will be assigned to the transmitting
% user
snr_meas=zeros(SP.N_frame,SP.N_SCFDMAsymb);
idx_symb_tot=0;
first_ref_symbol=1;

%-------------------------------------------------------------------------%
%                              Simulatore
%-------------------------------------------------------------------------%

clear PRM_che_results
Data=DataGen(SP);

SP.N_frame=1;

clear PRM_che
SNRcorr=mean(SP.SNR);
Eb_tot_sim=0;
No_tot_sim=0;
        
%-------------------------------------------------------------------------%
%                          Channel Preparation      
%-------------------------------------------------------------------------%
        
durata=SP.N_frame*SP.T_frame;
OR=SP.N_SCFDMAsymb/SP.T_frame;
       
        for iant_rx=1:SP.n_ant_RX
            switch SP.chan_type
                case 'AWGN'
                    % .' -> transpose matrix
                    paths=ones(SP.N_frame*SP.N_SCFDMAsymb,1).';
                    tau=0;

                case 'RAYL'
                    [tau,paths]=ext_channel('RAYL',NaN,SP.f0,SP.f_m,durata,OR);

                % 'EPA','EVA','ETU'
                case 'EPA'
                    [tau,paths]=ext_channel('EPA',NaN,SP.f0,SP.f_m,durata,OR);

                case 'EVA'
                    [tau,paths]=ext_channel('EVA',NaN,SP.f0,SP.f_m,durata,OR);

                case 'ETU'
                    [tau,paths]=ext_channel('ETU',NaN,SP.f0,SP.f_m,durata,OR);

                otherwise
                    error('Valore chan_type errato!')
            end
            
            tau=tau*1e-9; % [s]
            i_tau = ceil(SP.fs*tau*SP.f_upsamp)+1;
            
            for i1=1:length(i_tau)  
                h_chan(1:SP.N_frame*SP.N_SCFDMAsymb,i_tau(i1),iant_rx) = paths(i1,:).';
            end
            
            PRM_che.Lid = max(i_tau);
        end % iant_rx

        ichn=0;
        mpath_ext=[]; % interference from the previous symbol
        pattern_che=0;
        zero_len = SP.FFT_len - SP.Msymb;
        PRM_che.zero_offset=zero_len/2;
        
            idx_frame=1;

            frame_TX = [];
            frame_RX = [];
            frame_TX_chn = [];
            ifsgn = 0;
            k=1;
            j=1;            
            for idx_symb = 1:SP.N_SCFDMAsymb
                
                %-----------------------------------------------------------------%
                %                        Modulator SCFDMA
                %-----------------------------------------------------------------%
                
                % Note
                % The entire bandwidth is used for a single user
                
                % Generation of data bits
                DATA_TX(:,idx_symb) = Data(:,idx_symb);
                  
                % Mapping of bits in the costellation / reference symbols
                [SYMB_TX(:,idx_symb),reference_symbol,pilots_corr] = ConstMap(DATA_TX(:,idx_symb),SP,idx_frame,idx_symb);

                % M-DFT
                x_TX = fft(SYMB_TX(:,idx_symb).');
                
                if (reference_symbol)
                    if (first_ref_symbol)
                        SP.factor_ETX = round(sqrt(x_TX*x_TX'/length(x_TX))*100)/100;
                        first_ref_symbol=0;
                    end
                    x_TX = SP.factor_ETX * pilots_corr.';
                end
                
                % Subcarrier mapping
                x_TXmapped = [zeros(1,(zero_len)/2), x_TX, zeros(1,(zero_len)/2)];
                
                % N-IDFT
                y_TX = ifft(x_TXmapped);
                
                % Addition of the cyclic prefix
                y_TX_CP = AddCP(y_TX,idx_symb,SP);
                
                % Upsampling to simulate continuos time
                symbSCFDMA_TX = resample(y_TX_CP,SP.f_upsamp,1);
                symbTX_len(idx_symb) = length(symbSCFDMA_TX);
                
                % Vector of transmitted symbols
                frame_TX(1,ifsgn+1:ifsgn+length(symbSCFDMA_TX)) = symbSCFDMA_TX;
                
                %---------------------------------------------------------------------%
                %                               Channel
                %---------------------------------------------------------------------%
                
                % Channel index
                ichn=ichn+1;
                
                for iant_rx=1:SP.n_ant_RX
                    symbSCFDMA_TX_chn = conv(symbSCFDMA_TX, squeeze(h_chan(ichn,:,iant_rx)));
                    symbSCFDMA_TX_chn = symbSCFDMA_TX_chn(1:length(symbSCFDMA_TX)); 
                    % Multipath creation which modify the signal both in
                    % phase and amplitude
                    rho=0.2+(.2-.01)*rand(1,length(symbSCFDMA_TX_chn));
                    theta=0+(2*pi-0)*rand(1,length(symbSCFDMA_TX_chn));
                    mpath_ext=(rho.*exp(-2*pi*1i*theta));

                    if not(isempty(mpath_ext))
                          symbSCFDMA_TX_chn(1:length(symbSCFDMA_TX)) = ...
                            symbSCFDMA_TX_chn(1:length(symbSCFDMA_TX)) + mpath_ext(1:length(symbSCFDMA_TX));
                    end
                    
                    H_che_tmp=fft(padarray(squeeze(h_chan(ichn,iant_rx)),SP.FFT_len*SP.f_upsamp));
                    H_che_id(1:SP.Msymb,ichn,iant_rx)= ...
                       H_che_tmp(PRM_che.zero_offset+1:PRM_che.zero_offset+SP.Msymb).';
              
                    frame_TX_chn(iant_rx,ifsgn+1:ifsgn+length(symbSCFDMA_TX)) = symbSCFDMA_TX_chn;
                    
                end% iant_rx
                
                pattern_che=pattern_che+1;
                ifsgn = ifsgn + length(symbSCFDMA_TX);
                
            end % idx_symb

                        
            %---------------------------------------------------------------------%
            %                               SNR
            %---------------------------------------------------------------------%
            
            deltaIDFT = 10*log10(SP.Msymb/SP.FFT_len);
            deltaCP = 10*log10(length(y_TX)/length(y_TX_CP));
            snr = SNRcorr+deltaIDFT+deltaCP;
            
            for iant_rx=1:SP.n_ant_RX
                
                switch SP.SNR_comp
                    
                    case 'MEAS'
                        frame_RX(iant_rx,1:size(frame_TX_chn,2)) = ...
                            awgn(frame_TX_chn(iant_rx,:),snr,'measured');
                        
                    case '3GPP'
                        % SNR = S/N over a frame
                        S=frame_TX_chn(iant_rx,:)*frame_TX_chn(iant_rx,:)'/size(frame_TX_chn,2);
                        sigman=sqrt(S/(10^(snr/10)));
                        frame_RX(iant_rx,1:size(frame_TX_chn,2)) = ...
                            frame_TX_chn(iant_rx,1:size(frame_TX_chn,2)) + sqrt(2)/2*sigman*(randn(1,size(frame_TX_chn,2))+1i*randn(1,size(frame_TX_chn,2)));

                    otherwise
                        
                        frame_RX(iant_rx,1:size(frame_TX_chn,2)) = frame_TX_chn(iant_rx,1:size(frame_TX_chn,2));
                end
            end % iant_rx
            
            % Analysis of interfering frames for each used
            % type (number) of antennas  
            Frame_Interfering(nue_int).Ant1=frame_RX(1,:);
            Frame_Interfering(nue_int).Ant2=frame_RX(1:2,:);
            Frame_Interfering(nue_int).Ant4=frame_RX(1:4,:);
            Frame_Interfering(nue_int).Ant8=frame_RX(1:8,:);

            ue(nue_int).flag_tx_int=1;