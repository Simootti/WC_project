function log_flag = SimLog(SP,SER,BER,SNR_medio) 

%-------------------------------------------------------------------------%
%                              SimLog.m
%-------------------------------------------------------------------------%
%   Function which generates the log file of the simulation. 
% 
%   INPUT:
%    - SP             Simulator parameters (type: struct)
%    - SER            SER estimation
%    - BER            BER estimation
% 
%   OUTPUT:
%    - log_flag       0 = file created, -1 = file not created 
%    - logFile.txt    Log file of the simulation.
%     
%   Autore: Puccitelli Marco
%-------------------------------------------------------------------------%

    tmp_cd = cd;
    cd(SP.OutputDir);

    NtotSY = SP.Msymb*SP.N_SCFDMAsymb*SP.N_frame;
    NtotBIT = SP.mod_bitsymb*SP.Msymb*SP.N_SCFDMAsymb*SP.N_frame;

    fid = fopen('fileLog.txt','wt');

    fprintf(fid,'\n ====================================================== ');
    fprintf(fid,'\n              Simulatore SCFDMA - LogFile               ');
    fprintf(fid,'\n ====================================================== ');

    fprintf(fid,'\n\n\n Configurazione del simulatore: ');
    fprintf(fid,'\n ------------------------------ \n ');
    fprintf(fid,'\n Modulazione = %s \n', SP.modulazione);
    fprintf(fid,'\n Durata prefisso ciclico = %d [s] \n', SP.CP);
    fprintf(fid,'\n Durata simulazione = %d [s] \n', SP.T_sim);
    fprintf(fid,'\n Banda TX = %d [MHz] \n', SP.TX_bw);
    fprintf(fid,'\n Filtro trasmissione = %d  [1=rect/2=coseno rialzato] \n', SP.P_shape);
    fprintf(fid,'\n Fattore di upsampling = %d \n', SP.f_upsamp);


    fprintf(fid,'\n\n\n Configurazione del canale: ');
    fprintf(fid,'\n ------------------------------ \n ');
    switch SP.chan_type
        case 1
            fprintf(fid,'\n Canale = AWGN \n');
        case 2
            fprintf(fid,'\n Canale = Multipath + AWGN \n');
            fprintf(fid,'\n Ampiezza coeff. = %d \n', SP.chan_ampl);
            fprintf(fid,'\n Tau = %d \n', SP.chan_tau);
        otherwise
    end
    fprintf(fid,'\n SNR = %d [dB] \n', SP.SNR);

    fprintf(fid,'\n\n\n Output del simulatore: ');
    fprintf(fid,'\n ---------------------- \n');
    fprintf(fid,'\n #Simboli trasmessi = %d \n', NtotSY);
    fprintf(fid,'\n #Bit trasmessi = %d \n', NtotBIT);
    fprintf(fid,'\n SER = %d \n', SER);
    fprintf(fid,'\n BER = %d \n', BER);
    fprintf(fid,'\n SNR medio (misurato) = %d [dB] \n\n', SNR_medio);

    log_flag = fclose(fid);
    cd(tmp_cd);

%End_Of_Function
end

%End_Of_File