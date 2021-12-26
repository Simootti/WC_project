function [hfig,leg_string]=fig_LTE_uplink(file_save,type,hfig,str_flag,directory)
%
% [hfig,leg_string]=fig_LTE_uplink(file_save,type,hfig,directory)
%
% file_save
% type = 1 -> SER
%        2 -> BER

if (nargin==0)
    file_save='all';
    type=2;
    hfig=[];
    str_flag='';
    directory='.\results\';
end
if (nargin==1)
    type=2;
    hfig=[];
    str_flag='';
    directory='.\results\';
end
if (nargin==2)
    hfig=[];
    str_flag='';
    directory='.\results\';
end
if (nargin==3)
    str_flag='';
    directory='.\results\';
end
if (nargin==4)
    directory='.\results\';
end

leg_string=[];
if isempty(hfig)
    hfig=figure;
    grid on
    box on
    hold on
    new_figure=1;
else
    figure(hfig);
    new_figure=0;
end
%
if strcmp(file_save,'all')
    file_save=dir(directory);
end
if ischar(file_save)
    str_tmp=file_save;
    clear file_save
    if (length(str_tmp)>4)
        if not(strcmp(str_tmp(end-3:end),'.mat'))
            str_tmp=[str_tmp,'.mat'];
        end
    else
        str_tmp=[str_tmp,'.mat'];
    end
    file_save(1).name=str_tmp;
end

color_list=['k';'b';'r';'g'];
str_flag_list=['-+';'-*';'-o';'-s';'-d';'->';'-<';'-x'];
icl=1;
isf=0;

for ifl=1:length(file_save)
    load_file=0;
    if (length(file_save(ifl).name)>4)
        if (strcmp(file_save(ifl).name(end-3:end),'.mat'))
            load_file=1;
        end
    end
    if load_file&&(exist([directory,file_save(ifl).name],'file')==2)
        s=load([directory,file_save(ifl).name]);
        leg_string_corr=['SC-FDMA, ',s.SP.modulazione,', BW = ',num2str(s.SP.TX_bw),', SNR = ',s.SP.SNR_comp,', CH = ',s.SP.chan_type,', f_D = ',num2str(s.SP.f_m),', CHE = ',s.SP.CHE,', N_A(RX) = ',num2str(s.SP.n_ant_RX)];
        leg_string=strvcat(leg_string,leg_string_corr);

        if isempty(str_flag)
            isf=isf+1;
            if (isf==size(str_flag_list,1)+1)
                icl=icl+1;
                if (icl==size(color_list,1)+1)
                    icl=1;
                    isf=1;
                end
            end
            str_flag_corr=[color_list(icl,:),str_flag_list(isf,:)];
        else
            isf=isf+1;
            if (isf==size(str_flag,1)+1)
                isf=1;
            end
            str_flag_corr=[str_flag(isf,:)];
        end
        
        if (type==1)
            % Visualization of the SER graphic 
            semilogy(s.SP.SNR,mean(s.RES.SER,2)',str_flag_corr);

            %     ser_qpsk = erfc(sqrt(10.^((0:20)/10)/2)) ...
            %         - 1/4*erfc(sqrt(10.^((0:20)/10)/2)).^2 ;
            %     semilogy((0:20),ser_qpsk,'-');
            %     
            %     ser_16qam = 3/2*erfc(sqrt(10.^((0:28)/10)/10)) ...
            %         - 9/16*erfc(sqrt(10.^((0:28)/10)/10)).^2;
            %     semilogy((0:28),ser_16qam,'r-');
            %     
            %     ser_64qam = 7/4*erfc(sqrt(10.^((0:28)/10)/42))...
            %         - 49/64*erfc(sqrt(10.^((0:28)/10)/42)).^2;
            %     semilogy((0:28),ser_64qam,'g-');

            title(['SER: ',s.SP.modulazione]);
            xlabel('SNR [dB]');
            ylabel('P_S');
            % legend('SER misurato','QPSK','16-QAM','64-QAM'),
            ylim([1e-6 1e1]);

        elseif (type==2)

            semilogy(s.SP.SNR,mean(s.RES.BER,2)',str_flag_corr);

            %     ber_qpsk = 1/2*erfc(sqrt(10.^((0:20)/10)/2)) ...
            %         - 1/8*erfc(sqrt(10.^((0:20)/10)/2)).^2 ;
            %     semilogy((0:20),ber_qpsk,'-');
            %     
            %     ber_16qam = 3/8*erfc(sqrt(10.^((0:28)/10)/10)) ...
            %         - 9/64*erfc(sqrt(10.^((0:28)/10)/10)).^2;
            %     semilogy((0:28),ber_16qam,'r-');
            %     
            %     ber_64qam = 7/24*erfc(sqrt(10.^((0:28)/10)/42))...
            %         - 49/384*erfc(sqrt(10.^((0:28)/10)/42)).^2;
            %     semilogy((0:28),ber_64qam,'g-');

            title(['BER: ',s.SP.modulazione]);
            xlabel('SNR [dB]');
            ylabel('P_B');
            % legend('BER misurato','QPSK','16-QAM','64-QAM'),
            ylim([1e-6 1e0]);

        end

    end
    if (length(file_save)>1)&&(new_figure)
        legend(leg_string);
    end
end % ifl
