%% Function to print the graphs

EbN_dB = SP.SNR;
EbNo=10.^(0.1*EbN_dB);
NR = n_ant_RX_vect; % number of receiver antennas

switch combiner
    
    case 'MRC'

        % practical MRC combiner
        x = input("In order to plot the BER with N user interfering, digit number N, press enter to finish\n");
        q=1;
        L={};
        
        while not(isempty(x))
            k = input("Enter number of receiver antennas, press enter to stop\n");
            while not(isempty(k))
                txt = [num2str(k),'-ant MRC of user ',num2str(x)];
                switch k 
                    case 1
                        semilogy (EbN_dB,BER_PER_GRAFICO_CON_INTERFERENZA(x+1).ber1,'linewidth',1.0)
                    case 2
                        semilogy (EbN_dB,BER_PER_GRAFICO_CON_INTERFERENZA(x+1).ber2,'linewidth',1.0)
                    case 4
                        semilogy (EbN_dB,BER_PER_GRAFICO_CON_INTERFERENZA(x+1).ber4,'linewidth',1.0)
                    case 8
                        semilogy (EbN_dB,BER_PER_GRAFICO_CON_INTERFERENZA(x+1).ber8,'linewidth',1.0)
                end
                hold on
                L{q}=[num2str(k),'-ant MRC with ',num2str(x),' users int.'];
                q=q+1;
                
                k = input("Enter number of receiver antennas, press enter to stop\n");
            end
             x = input("In order to plot another BER with N user interfering, digit number N, press enter to finish\n");  
        end
        legend(L);
        title({'QPSK BER di un MP fading channel','[1,2,4,8] ANTENNAS con MRC combiner'});       
        axis([0 15 10^-5 1]);
        grid on
        xlabel('EbNo(dB)')
        ylabel('BER')
  
    case 'MAX_SEL'

        % practical MAX_SEL combiner
        x = input("In order to plot the BER with N user interfering, digit number N, press enter to finish\n");
        q=1;
        L={};
        
        while not(isempty(x))
            k = input("Enter number of receiver antennas, press enter to stop\n");
            while not(isempty(k))
                txt = [num2str(k),'-ant MS of user ',num2str(x)];
                switch k 
                    case 1
                        semilogy (EbN_dB,BER_PER_GRAFICO_CON_INTERFERENZA(x+1).ber1,'linewidth',1.0)
                    case 2
                        semilogy (EbN_dB,BER_PER_GRAFICO_CON_INTERFERENZA(x+1).ber2,'linewidth',1.0)
                    case 4
                        semilogy (EbN_dB,BER_PER_GRAFICO_CON_INTERFERENZA(x+1).ber4,'linewidth',1.0)
                    case 8
                        semilogy (EbN_dB,BER_PER_GRAFICO_CON_INTERFERENZA(x+1).ber8,'linewidth',1.0)
                end
                hold on
                L{q}=[num2str(k),'-ant MS with ',num2str(x),' users int.'];
                q=q+1;
                
                k = input("Enter number of receiver antennas, press enter to stop\n");
            end
             x = input("In order to plot another BER with N user interfering, digit number N, press enter to finish\n");  
        end
        legend(L);   
        title({'QPSK BER di un MP fading channel','[1,2,4,8] ANTENNAS con MAX SELECTION combiner'});
        axis([0 15 10^-5 1]);
        grid on
        xlabel('EbNo(dB)')
        ylabel('BER')
end