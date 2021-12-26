
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     SCFDMA Simulator - Uplink LTE
%%                    Base Station parameters ln 177                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Simulator in base band of an LTE uplink system.
%   The transmettitor generates the signal for a period equal to the 
%   subframe duration (1ms) and executes the following operations:
%       1) Generation of data bits
%       2) Mapping of bits in the chosen costellation (QPSK,16QAM,64QAM)
%       3) Execution of the DFT and subcarrier mapping (type: "localized")
%       4) Execution of the IDFT 
%       5) Addition of the cyclic prefix
%       6) Upsampling to simuluate the continuos time 
%   
%   The subframe transmission happens on a channel definied by the user:
%       1) AWGN channel 
%       2) "Multipath + AWGN" channel
%       NB: The max latency that can be simulated is equal to: t_frame = 1e-3 [s]
% 
%   The receiver executes the following operations:
%       1) Downsampling to simuluate the sampling
%       2) Removal of the cyclic prefix
%       3) Execution of the DFT and subcarrier demapping 
%       4) Execution of the IDFT
%       5) Decisor of "Hard-decision" type on received symbols
%       6) Demapping of the costellation to obtain the transmitted bit
% 
%   At each iteration, the transmitter evaluate the PAPR (Peak to 
%   average power ratio) and the receiver evaluate the SNR
%   measured for each transmitted SCFDMA symbol
%   and the number of wrong bit/symbols. 
%   Finally curves of BER, SER, measured SNR and received symbols 
%   (if in the cycle the "plot flag" is placed equal to 1) are visualized.
% 
%   The logfile and the matrix of the parameters of the simulation are
%   automatically saved in a folder inside the current process,
%   while in order to save input and output data is needed to activate
%   the option (simP.Sv_flag). 
% 
%   References: 3GPP TS 36.211 V8.8.0
%               3GPP TS 36.101 V9.5.0   
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
%                          Parametri simulatore 
%-------------------------------------------------------------------------%

% Parametri Modulatore/Demodulatore
if not(exist('Mod_type'))
    simP.Mod_type = 1;               % 1=QPSK , 2=16QAM , 3=64QAM   
else
    simP.Mod_type = Mod_type;        % 1=QPSK , 2=16QAM , 3=64QAM   
end
simP.CP_type = 1;                    % Cyclic Prefix: 1=Normal o 2=Extended
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
    simP.LSi_red_factor = 20;        % 0, 1 (ideal L)
else
    simP.LSi_red_factor = LSi_red_factor;
end

simP.chan_ampl = [1 .2];             % Channel coefficients vector (2 rays)
simP.chan_tau = [0 5e-6];            % Vector of channel latencies (tau) [s]
simP.f0 = 2e+9;                      % carrier frequency, 2GHz
if not(exist('f_m'))
    simP.f_m = 5;                    % Max Doppler
else
    simP.f_m = f_m;
end

%% Parametri BS
if not(exist('n_ant_RX'))
     simP.n_ant_RX = 1;
else
    simP.n_ant_RX = n_ant_RX;
end

% Simulation parameters
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
%                              Simulator
%-------------------------------------------------------------------------%

h = waitbar(0,'Attendere...');
tic

clear PRM_che_results
Data=DataGen(SP);

SP.N_frame=M_frame;

for isnr=1:length(SP.SNR)

    clear PRM_che
    SNRcorr=SP.SNR(isnr);    
    for isim=1:SP.N_sim
            Eb_tot_sim=0;
            No_tot_sim=0;
        
%-------------------------------------------------------------------------%
%                         Channel preparation              
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
        
        for idx_frame = 1:SP.N_frame

            frame_TX = [];
            frame_RX = [];
            frame_TX_chn = [];
            ifsgn = 0;
            k=1;
            j=1;            
            for idx_symb = 1:SP.N_SCFDMAsymb
                
                %-----------------------------------------------------------------%
                %                        SCFDMA Modulator
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
            
            % Measure of Peak-to-Average Power Ratio
            RES.PAPR(isnr,isim)=RES.PAPR(isnr,isim)+(max(abs(frame_TX).^2)/mean(abs(frame_TX).^2));
                        
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
                        %error('No noise set');
                end
            end % iant_rx
            frame_tot_interfering=zeros(size(frame_RX,1),size(frame_RX,2));     
            
            
            %-----------------------------------------------------------------%
            %                       Receiver SCFDMA
            %-----------------------------------------------------------------%
            %-----------------------------------------------------------------%
            %                       Channel estimation
            % Notes:
            % CHE is performed independently for each frame 
            % (To do : extend to generale pattern = SP.Pilots)
            %-----------------------------------------------------------------%
            
            for iant_rx=1:SP.n_ant_RX
                % Buffer generation
                Buffer=zeros(SP.Msymb,SP.N_SCFDMAsymb);
                Buffern=zeros(SP.Msymb,SP.N_SCFDMAsymb); % no noise
                Diff_cp_RX=[];
                Unused_band_noise=0;
                
                for idx_symb = 1:SP.N_SCFDMAsymb
                    % Received symbol
                    if idx_symb == 1
                        start_sample = 1;
                        end_sample = start_sample+symbTX_len(idx_symb)-1;
                    else
                        start_sample = 1 + sum(symbTX_len(1:idx_symb-1));
                        end_sample = start_sample+symbTX_len(idx_symb)-1;
                    end           
                    
                    symbSCFDMA_RX = frame_RX(iant_rx,start_sample:end_sample);
                    
                    % Downsampling in order to simulate the sampling
                    y_RX_CP = symbSCFDMA_RX(1:SP.f_upsamp:end);
                    
                    % Removal of the cyclic prefix
                    [y_RX, cp_RX, len_cp_RX] = RemCP(y_RX_CP,idx_symb,SP);
                    
                    % N-DFT
                    x_RXmapped = fft(y_RX);

                    Buffer(:,idx_symb) = x_RXmapped(zero_len/2+1:end-(zero_len/2)).';
                    Unused_band_noise = ...
                        Unused_band_noise +...
                           mean(abs(x_RXmapped([1:zero_len/2,length(x_RXmapped)-(zero_len/2)+1:length(x_RXmapped)])).^2);
                    
                    if (SP.noise_comp==1)
                        symbSCFDMA_RXn = frame_TX_chn(iant_rx,start_sample:end_sample);
                        y_RX_CPn = symbSCFDMA_RXn(1:SP.f_upsamp:end);
                        y_RXn = RemCP(y_RX_CPn,idx_symb,SP);
                        x_RXmappedn = fft(y_RXn);
                        Buffern(:,idx_symb) = x_RXmappedn(zero_len/2+1:end-(zero_len/2)).';
                    end
                    
                end % idx_symb

                Diff_cp_RX = abs(Diff_cp_RX/SP.N_SCFDMAsymb).^2;
                Unused_band_noise = Unused_band_noise/SP.N_SCFDMAsymb;
                pattern_che=0;
                
                % Estimation of NV + L
                switch SP.NVLE
                    
                   case 'IDEAL'
                        PRM_che.sigma2nf(iant_rx,1:size(Buffer,1))=mean(abs(Buffer-Buffern).^2,2);
                        PRM_che.sigma2n(iant_rx)=mean(PRM_che.sigma2nf(iant_rx,:));
                        PRM_che.L(iant_rx)=PRM_che.Lid;
                    
                    case 'CP'
                        % Analysis of Diff_cp_RX
                        sigma2n_est=2*Unused_band_noise;
                        PRM_che.sigma2n(iant_rx)=sigma2n_est;
                        Pfa=0.1;
                        zthr=-sigma2n_est*log(1-(1-Pfa));
                        L0=find(Diff_cp_RX>zthr);
                        if isempty(L0)
                            PRM_che.L(iant_rx)=1;
                        else
                            PRM_che.L(iant_rx)=1+max(L0);
                        end
                end
                
                % Channel estimation
                % Ideal reference
                a=1;
                if (idx_frame>1)
                    [H_che_id(:,ichn+1:ichn+size(CHE_ofdm('LSf',SP,Buffern,PRM_che),2),iant_rx)] = ... 
                        CHE_ofdm('LSf',SP,Buffern,PRM_che);
                else
                    [H_che_id(:,1:ichn,iant_rx)] = CHE_ofdm('LSf',SP,Buffern,PRM_che);
                end
                
                %Switch of the channel estimation set in Main
                switch SP.CHE
                    case 'LS'
                        [H_che(:,:,iant_rx),PRM_che] = CHE_ofdm('LS',SP,Buffer,PRM_che);
                        
                    case 'LSf'
                        [H_che(:,:,iant_rx),PRM_che] = CHE_ofdm('LSf',SP,Buffer,PRM_che);
                        
                    case 'MMSE'
                        [H_che(:,:,iant_rx),PRM_che] = CHE_ofdm('MMSE',SP,Buffer,PRM_che);
                        
                    case 'LSi'
                        [H_che(:,:,iant_rx),PRM_che] = CHE_ofdm('LSi',SP,Buffer,PRM_che);
                        
                    case 'LSb'
                        [H_che(:,:,iant_rx),PRM_che] = CHE_ofdm('LSb',SP,Buffer,PRM_che);
                        
                    case 'IDEAL'
                        H_che(:,:,iant_rx)=H_che_id(:,:,iant_rx);
                        
                    otherwise
                        H_che=ones(size(Buffer));
                        warning('No CHE performed');
                end
            end % iant_rx
   
            %% Interference set in the Main
            if(INTERFERENCE)
                if ue(nue).number_users_interfering
                    for k=1:numel(ue(nue).u_int)
                         idx_int=ue(nue).u_int(k);
                         switch n_ant_RX
                             case 1
                                 frm_int=Frame_Interfering(idx_int).Ant1;
                             case 2
                                 frm_int=Frame_Interfering(idx_int).Ant2;
                             case 4
                                 frm_int=Frame_Interfering(idx_int).Ant4;
                             case 8
                                 frm_int=Frame_Interfering(idx_int).Ant8;
                         end

                         weight_interference=sqrt(dBm_to_watt(Bs(ue(nue).NNCellID).Pr_nominale(idx_int)))/...
                                         sqrt(dBm_to_watt(Bs(ue(nue).NNCellID).Pr_nominale(nue)));

                         frame_tot_interfering = frame_tot_interfering +...
                             (frm_int).*weight_interference;
                    end
                end
            end

            % Addition of the interference
            frame_RX=frame_RX+frame_tot_interfering;

            % For MAX SELECTION COMBINER EQUALIZATION is needed, so before
            % combining we consider the best channel gains, i.e. those with
            % with the highest power
            Channel_Gains_sum=0;
            Gains=0;
        
            %% ------------------------------------------------------------%
            %                    Demodulator SCFDMA                       %
            %%-------------------------------------------------------------%

            for idx_symb = 1:SP.N_SCFDMAsymb
                
                % Received symbol
                if idx_symb == 1
                    start_sample = 1;
                    end_sample = start_sample+symbTX_len(idx_symb)-1;
                else
                    start_sample = 1 + sum(symbTX_len(1:idx_symb-1));
                    end_sample = start_sample+symbTX_len(idx_symb)-1;
                end

                x_RX = zeros(1,SP.Msymb);
                x_RX_nn = zeros(1,SP.Msymb);
                N_RX = zeros(1,SP.Msymb);
                Signal_antenna={};
                power_signal=0;
                
                for iant_rx=1:SP.n_ant_RX
                    symbSCFDMA_RX = frame_RX(iant_rx,start_sample:end_sample);
                    symbSCFDMA_RX_nn = frame_TX_chn(iant_rx,start_sample:end_sample);
                    
                    % Downsampling to simulate the sampling
                    y_RX_CP = symbSCFDMA_RX(1:SP.f_upsamp:end);
                    y_RX_CP_nn = symbSCFDMA_RX_nn(1:SP.f_upsamp:end);
                
                    % Removal of the cyclic prefix
                    y_RX = RemCP(y_RX_CP,idx_symb,SP);
                    y_RX_nn = RemCP(y_RX_CP_nn,idx_symb,SP);
                
                    % N-DFT
                    x_RXmapped = fft(y_RX);
                    x_RXmapped_nn = fft(y_RX_nn);
                
                    % Subcarrier demapping
                    x_RX_ant = x_RXmapped(zero_len/2+1:end-(zero_len/2));
                    x_RX_ant_nn = x_RXmapped_nn(zero_len/2+1:end-(zero_len/2));

                    Signal_antenna(iant_rx).sgnl=x_RX_ant;
                    Signal_antenna(iant_rx).sgnl_nn=x_RX_ant_nn;
                    
                    %% ------------------------------------------------------------%
                    %                          Combiner                            %
                    %%-------------------------------------------------------------%
                    
                     switch combiner   
                         
                         case 'MRC'
                             
                             % Channel Estimation for each antenna
                             % Hcorr_2 has to be conjugated (apex ')
                             H_corr2 = squeeze(H_che(:,idx_symb,iant_rx))';
                            
                             %Hcorr_2 has been conjugated before
                             x_RX = x_RX + x_RX_ant.*(H_corr2);
                             
                             % As the definition of MRC combiner, the demapped 
                             % signal can be combined with channel estimation,
                             % each signal is by the gains of the corrensponding
                             % antenna
                             x_RX_nn = x_RX_nn + x_RX_ant_nn.*(H_corr2);
                             N_RX = N_RX + abs(H_corr2).^2;
                             
                         case 'MAX_SEL'
                             
                             H_corr2 = squeeze(H_che(:,idx_symb,iant_rx))';
                             if power_signal<sqrt(sum(abs(H_corr2).^2))
                                 power_signal=sqrt(sum(abs(H_corr2).^2));
                                 signal=x_RX_ant;
                                 signal_nn=x_RX_ant_nn;   
                                 gains=H_corr2;
                             end
                     end
                end % iant_rx
                
                % Equalization needed for max_sel, but in this
                % case we need to use the 'best' channel gains,
                % i.e. the most powerful signal
                switch combiner
                    case 'MAX_SEL'
                        x_RX = signal.*gains;
                        x_RX_nn = signal_nn.*gains;
                        N_RX = abs(gains).^2;                    
                end
 
                x_RX = x_RX ./ N_RX;
                x_RX_nn = x_RX_nn ./ N_RX;
                
                % M-IDFT
                SYMB_RX(:,idx_symb) = ifft(x_RX);
                SYMB_RX_nn(:,idx_symb) = ifft(x_RX_nn);

                % Decisore
                SYMB_stima(:,idx_symb) = Decisor(SYMB_RX(:,idx_symb).',SP);

                % Demapping of the symbol's costellation
                [DATA_RX(:,idx_symb), TYPE_DATA, TYPE_SYM] =...
                    ConstDemap(SYMB_stima(:,idx_symb).',SP,idx_frame,idx_symb);

                %-----------------------------------------------------------------%
                %           Measure and evaluation of wrong symbols/bit 
                %-----------------------------------------------------------------%                 
                
                % Comparison Symbols and Bits transmitted/received
                i_bit = find(TYPE_DATA==1);
                Nbit = length(i_bit);
                
                if (Nbit>0)
                    i_sy = find(TYPE_SYM==1);
                    Nsy = length(i_sy); % round(Nbit/SP.mod_bitsymb); 
                    [Nsy_err, S_ratio, Sy_err] = symerr(SYMB_TX(i_sy,idx_symb),SYMB_stima(i_sy,idx_symb));
                    [Nbit_err, B_ratio, Bit_err] = biterr(DATA_TX(i_bit,idx_symb),DATA_RX(i_bit,idx_symb));

                    % Updating the "error counters"
                    RES.SEcount(isnr,isim) = RES.SEcount(isnr,isim) + Nsy_err;
                    RES.BEcount(isnr,isim) = RES.BEcount(isnr,isim) + Nbit_err;
                    RES.Scount(isnr,isim) = RES.Scount(isnr,isim) + Nsy;
                    RES.Bcount(isnr,isim) = RES.Bcount(isnr,isim) + Nbit;
                end
                
                const_plot=0;
                if const_plot
                    % Visualization of the costellation of received symbols
                    Output.SY_RX = SYMB_RX(:,idx_symb);
                    Output.SY_err = Sy_err;
                    constIQplot(SP,Output);
                    hold off,
                end
                
                if (SP.Sv_flag == 2)&&(isnr==1)&&(isim==1)
                    % Savings of data and symbols
                    Input.DATA = [Input.DATA DATA_TX];
                    Input.SY_TX = [Input.SY_TX SYMB_TX];
                    Output.SY_RX = [Output.SY_RX SYMB_RX];
                    Output.SY_stima = [Output.SY_stima SYMB_stima];
                    Output.SY_err = [Output.SY_err Sy_err];
                    Output.DATA = [Output.DATA DATA_RX];
                end
                
                idx_symb_tot=idx_symb_tot+1;           

                % Updating the "state bar"
                waitbar(idx_symb_tot/(SP.N_SCFDMAsymb*SP.N_sim*length(SP.SNR)*SP.N_frame),h,['SNR = ',...
                    num2str(SNRcorr),' dB - sim (',num2str(isim),'/',...
                        num2str(SP.N_sim),') - subframe (' ...
                            num2str(idx_frame),'/',num2str(SP.N_frame),') - ',num2str(n_ant_RX),' antennas ']);
                                % Uncomment if we want to add the percentage on the waitbar
                                % num2str(ceil((idx_symb_tot/(SP.N_SCFDMAsymb*SP.N_sim*length(SP.SNR)*SP.N_frame))*100)) ' %']);
                
            end % idx_symb
            
            
        end % idx_frame
        
    end % isim

    % Store CHE results
    if isfield(PRM_che,'Len_filter_GB')
        PRM_che_results.Len_filter_GB(isnr)=PRM_che.Len_filter_GB;
    end         
end            

close(h);

%-------------------------------------------------------------------------%
%                           Performance Analysis
%-------------------------------------------------------------------------%

% Symbol Error Rate
RES.NtotSY = SP.Msymb*SP.N_SCFDMAsymb*SP.N_frame;
RES.SER = RES.SEcount./RES.Scount;

% Bit Error Rate
RES.NtotBIT = SP.mod_bitsymb*SP.Msymb*SP.N_SCFDMAsymb*SP.N_frame;
RES.BER = RES.BEcount./RES.Bcount;

RES.PAPR=10*log10(RES.PAPR/SP.N_frame);
RES.Err_che=RES.Err_che/SP.N_sim;

RES.deltaIDFT = deltaIDFT;
RES.deltaCP = deltaCP;

SP.plot_fig=0;
if (SP.plot_fig)
    % Visualization of the SER graph
    figure('Name','SER');
    semilogy(SP.SNR,mean(RES.SER,2)','k*'), grid, hold on,
    
    ser_qpsk = erfc(sqrt(10.^((0:15)/10)/2)) ...
        - 1/4*erfc(sqrt(10.^((0:15)/10)/2)).^2 ;
    semilogy((0:15),ser_qpsk,'-');
    
    ser_16qam = 3/2*erfc(sqrt(10.^((0:15)/10)/10)) ...
        - 9/16*erfc(sqrt(10.^((0:15)/10)/10)).^2;
    semilogy((0:15),ser_16qam,'r-');
    
    ser_64qam = 7/4*erfc(sqrt(10.^((0:15)/10)/42))...
        - 49/64*erfc(sqrt(10.^((0:15)/10)/42)).^2;
    semilogy((0:15),ser_64qam,'g-');
    
    title(['SER: ' SP.modulazione ', SNR = ' num2str(SP.SNR) 'dB']), 
    xlabel('SNR_d_B'), ylabel('P( \epsilon )'),
    legend('SER misurato','QPSK','16-QAM','64-QAM'),
    ylim([1e-6 1e1]);

    % Visualizza il grafico per BER
    figure('Name','BER'),
    semilogy(SP.SNR,mean(RES.BER,2)','k*'), grid, hold on,
    
    ber_qpsk = 1/2*erfc(sqrt(10.^((0:15)/10)/2)) ...
        - 1/8*erfc(sqrt(10.^((0:15)/10)/2)).^2 ;
    semilogy((0:15),ber_qpsk,'-');
    
    ber_16qam = 3/8*erfc(sqrt(10.^((0:15)/10)/10)) ...
        - 9/64*erfc(sqrt(10.^((0:15)/10)/10)).^2;
    semilogy((0:15),ber_16qam,'r-');
    
    ber_64qam = 7/24*erfc(sqrt(10.^((0:15)/10)/42))...
        - 49/384*erfc(sqrt(10.^((0:15)/10)/42)).^2;
    semilogy((0:15),ber_64qam,'g-');
    
    title(['BER: ' SP.modulazione ', SNR = ' num2str(SP.SNR) 'dB']),
    xlabel('SNR_d_B'), ylabel('P( \epsilon )'),
    legend('BER misurato','QPSK','16-QAM','64-QAM'),
    ylim([1e-6 1e1]);
    
end

% Saving the input/output file
if (SP.Sv_flag == 1)
    % Salva gli input
    file_input = [SP.OutputDir,'In_',SP.file_save,date,'.mat'];
    %save(file_input,'Input SP','-v6');
    % Save the output
    file_output = [SP.OutputDir,SP.file_save,'_',date,'.mat'];
    % The following road has been commented out to avoid errors
    % eval(['save ',file_output,' SP RES PRM_che_results -v6']);
end

if (SP.Sv_flag == 2)
    file_output = [SP.OutputDir,SP.file_save,'_',date,'.mat'];
    %eval(['save ',file_output,' SP RES PRM_che_results Input Output -v6']);
end

disp(' ')
disp('Simulazione terminata')
disp(' ')
toc

disp(' ')
disp(['Logfile salvato nella cartella: ' SP.OutputDir])
disp(' ')

if (SP.Sv_flag == 1)||(SP.Sv_flag == 2)    
    disp(['Dati salvati nella cartella: ' SP.OutputDir])
    disp(' ')
end

% End_Of_File