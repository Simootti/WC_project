function [H_CHE,PRM_CHE_new] = CHE_ofdm(mode,SP,Buffer,PRM_CHE)
%
%% CHE=Channel Estimation
% [H_CHE,PRM_CHE_new] = CHE_ofdm(mode,SP,Buffer)
%
% mode = 'LS'  - Least squares + linear interpolation
%      = 'LSf' - Frequency domain LS + linear interpolation
%      = 'LSi' - Approx frequency domain LS through interpolation + linear interpolation
%      = 'LSb' - Block Least squares + linear interpolation
%
% SP        = see SimConfig.m
% Buffer    = Frequency samples (same size of SP.Pilots)
% PRM_CHE   = Additional parameters for CHE
%             PRM_CHE can be updated by the function
%             .L = estimated length of the impulse response

%         SP.Pilots=zeros(SP.Msymb,SP.N_SCFDMAsymb);
%         SP.Pilot_symbol

% compute max length of cyclic prefix
switch SP.CP_type       
    case 1  % Normal CP 
        max_len_cp=round((SP.CP(2)*SP.fs));
    case 2  % Extended CP
        max_len_cp=round((SP.CP(1)*SP.fs));
end

H_CHE=ones(size(Buffer)); %[1 1 1 ... # samples in freq]
PRM_CHE_new=PRM_CHE;
if isempty(PRM_CHE_new)
    clear PRM_CHE_new
    PRM_CHE_new.L=[];
end
if not(isfield(PRM_CHE_new,'L'))||isempty(PRM_CHE_new.L)
    PRM_CHE_new.L=max_len_cp;
end
if not(isfield(PRM_CHE_new,'alpha'))||isempty(PRM_CHE_new.alpha)
    PRM_CHE_new.alpha=0.1;
end
if (SP.Version_CHE==1)
    L=PRM_CHE_new.Lid;
elseif (SP.Version_CHE==0)
    L=max_len_cp;
elseif (SP.Version_CHE==2)
    L=PRM_CHE_new.L;
end
pattern_len=size(SP.Pilots,2);

if (size(Buffer)==size(SP.Pilots))
    linear_interpolation=0;
    if strcmp(mode,'LS')
        %
        % Frequency estimation 
        li_taps=[];
        for isym=1:size(SP.Pilots,2)
            flag_pilots_corr=SP.Flag_Pilots(:,isym);
            i_pilots=find(flag_pilots_corr');
            if not(isempty(i_pilots))
                li_taps=[li_taps,isym];
                % h_LS on the current symbol
                N_pilots=length(i_pilots);
                build_GB=0;
                if not(isfield(PRM_CHE_new,'GB'))
                    build_GB=1;
                else
                    if not(isequal(size(PRM_CHE_new.BI,1),L))
                        build_GB=1;
                    end
                end
                if build_GB
                    B=zeros(N_pilots,L);
                    G=zeros(SP.FFT_len,L);
                    for iB=1:size(B,2)
                        B(:,iB)=exp(-2*1i*pi*(i_pilots-1+PRM_CHE_new.zero_offset)*(iB-1)/SP.FFT_len).';
                        % B(:,iB)=exp(-2*1i*pi*(i_pilots-1+0)*(iB-1)/SP.FFT_len).';
                        G(:,iB)=exp(-2*1i*pi*[0:SP.FFT_len-1]*(iB-1)/SP.FFT_len).';
                    end
                    PRM_CHE_new.BI=inv(B'*B+PRM_CHE_new.alpha*eye(L))*B'; % L x N_pilots
                    PRM_CHE_new.GB=G*PRM_CHE_new.BI;
                    PRM_CHE_new.GB=PRM_CHE_new.GB(PRM_CHE_new.zero_offset+1:PRM_CHE_new.zero_offset+size(Buffer,1),:);
                end
                Yp=Buffer(i_pilots,isym).*conj(SP.Pilots(i_pilots,isym))./(SP.factor_ETX*abs(SP.Pilots(i_pilots,isym)).^2);
                % h_LS=BI*Yp;
                % H_corr=fft(zero_pad(h_LS.',0,SP.FFT_len));
                % H_CHE(:,isym)=H_corr(PRM_CHE_new.zero_offset+1:PRM_CHE_new.zero_offset+size(Buffer,1)).';
                H_CHE(:,isym)=PRM_CHE_new.GB*Yp;
            end
        end
        linear_interpolation=1;
        %
    elseif strcmp(mode,'LSi')
        %
        % Frequency estimation 
        li_taps=[];
        for isym=1:size(SP.Pilots,2)
            flag_pilots_corr=SP.Flag_Pilots(:,isym);
            i_pilots=find(flag_pilots_corr');
            if not(isempty(i_pilots))
                li_taps=[li_taps,isym];
                % h_LS on the current symbol
                N_pilots=length(i_pilots);
                build_GB=0;
                if not(isfield(PRM_CHE_new,'filter_GB'))
                    build_GB=1;
                else
                    if not(isequal(size(PRM_CHE_new.BI,1),L))
                        build_GB=1;
                    end
                end
                if build_GB
                    B=zeros(N_pilots,L);
                    G=zeros(SP.FFT_len,L);
                    for iB=1:size(B,2)
                        B(:,iB)=exp(-2*1i*pi*(i_pilots-1+PRM_CHE_new.zero_offset)*(iB-1)/SP.FFT_len).';
                        G(:,iB)=exp(-2*1i*pi*[0:SP.FFT_len-1]*(iB-1)/SP.FFT_len).';
                    end
                    PRM_CHE_new.BI=inv(B'*B+PRM_CHE_new.alpha*eye(L))*B'; % L x N_pilots
                    PRM_CHE_new.GB=G*PRM_CHE_new.BI;
                    PRM_CHE_new.GB=PRM_CHE_new.GB(PRM_CHE_new.zero_offset+1:PRM_CHE_new.zero_offset+size(Buffer,1),:);
                    PRM_CHE_new.filter_GB=fliplr(PRM_CHE_new.GB(round(size(Buffer,1)/2),:));
                    PRM_CHE_new.ind0_filter_GB=length(PRM_CHE_new.filter_GB)-round(size(Buffer,1)/2);
                    [Htmp]=redh([PRM_CHE_new.ind0_filter_GB,PRM_CHE_new.filter_GB],SP.LSi_red_factor);
                    PRM_CHE_new.filter_GB=Htmp(1,2:end);
                    PRM_CHE_new.ind0_filter_GB=Htmp(1);
                    Href=conv(ones(1,N_pilots),PRM_CHE_new.filter_GB);
                    Href=Href(1,PRM_CHE_new.ind0_filter_GB:PRM_CHE_new.ind0_filter_GB+size(Buffer,1)-1).';
                    PRM_CHE_new.Kref=1./abs(Href);
                    PRM_CHE_new.Len_filter_GB=length(Htmp)-1;
                end
                Yp=Buffer(i_pilots,isym).*conj(SP.Pilots(i_pilots,isym))./(SP.factor_ETX*abs(SP.Pilots(i_pilots,isym)).^2);
                H_corr=conv(Yp.',PRM_CHE_new.filter_GB);
                H_CHE(:,isym)=PRM_CHE_new.Kref.*H_corr(1,PRM_CHE_new.ind0_filter_GB:PRM_CHE_new.ind0_filter_GB+size(Buffer,1)-1).';
            end
        end
        linear_interpolation=1;
        %
    elseif strcmp(mode,'LSb')
        %
        % Frequency estimate 
        li_taps=[];
        for isym=1:size(SP.Pilots,2)
            flag_pilots_corr=SP.Flag_Pilots(:,isym);
            i_pilots=find(flag_pilots_corr');
            if not(isempty(i_pilots))
                li_taps=[li_taps,isym];
                % h_LS on the current symbol
                N_pilots=length(i_pilots);
                build_GB=0;
                if not(isfield(PRM_CHE_new,'GB'))
                    build_GB=1;
                else
                    if not(isequal(size(PRM_CHE_new.BI,1),min(L,12)))
                        build_GB=1;
                    end
                end
                if build_GB
                    B=zeros(N_pilots,min(L,12));
                    G=zeros(SP.FFT_len,min(L,12));
                    GB=zeros(N_pilots,N_pilots);
                    for iB=1:size(B,2)
                        G(:,iB)=exp(-2*1i*pi*[0:SP.FFT_len-1]*(iB-1)/SP.FFT_len).';
                    end
                    nRB=N_pilots/12;
                    for inRB=1:nRB
                        indcorr=[(inRB-1)*12+1:inRB*12];
                        Bcorr=zeros(12,min(L,12));
                        for iB=1:size(Bcorr,2)
                            Bcorr(:,iB)=exp(-2*1i*pi*(i_pilots(indcorr)-1+PRM_CHE_new.zero_offset)*(iB-1)/SP.FFT_len).';
                            % B(:,iB)=exp(-2*1i*pi*(i_pilots-1+0)*(iB-1)/SP.FFT_len).';
                        end
                        BIcorr=inv(Bcorr'*Bcorr+PRM_CHE_new.alpha*eye(min(L,12)))*Bcorr'; % L x N_pilots
                        GBcorr=G*BIcorr;
                        GBcorr=GBcorr(PRM_CHE_new.zero_offset+i_pilots(indcorr),:); % N_pilots x N_pilots
                        GB(indcorr,indcorr)=GBcorr;
                    end
                    PRM_CHE_new.GB=GB;
                    PRM_CHE_new.BI=NaN*ones(min(L,12),N_pilots); % L x N_pilots
                end
                Yp=Buffer(i_pilots,isym).*conj(SP.Pilots(i_pilots,isym))./(SP.factor_ETX*abs(SP.Pilots(i_pilots,isym)).^2);
                % h_LS=BI*Yp;
                % H_corr=fft(zero_pad(h_LS.',0,SP.FFT_len));
                % H_CHE(:,isym)=H_corr(PRM_CHE_new.zero_offset+1:PRM_CHE_new.zero_offset+size(Buffer,1)).';
                H_CHE(:,isym)=PRM_CHE_new.GB*Yp;
            end
        end
        linear_interpolation=1;
        %
    elseif strcmp(mode,'LSf')
        %
        % Frequency estimate 
        li_taps=[];
        for isym=1:size(SP.Pilots,2)
            flag_pilots_corr=SP.Flag_Pilots(:,isym);
            i_pilots=find(flag_pilots_corr');
            if not(isempty(i_pilots))
                li_taps=[li_taps,isym];
                % h_LS on the current symbol
                N_pilots=length(i_pilots);
                Yp=Buffer(i_pilots,isym).*conj(SP.Pilots(i_pilots,isym))./(SP.factor_ETX*abs(SP.Pilots(i_pilots,isym)).^2);
                Hcorr=Yp;
                H_CHE(:,isym)=Hcorr.';
            end
        end
        linear_interpolation=1;
        %
    elseif strcmp(mode,'MMSE')
        disp('Not implemented yet');
        % ...
    else
        warning('Incorrect CHE mode');
    end
    
    % Linear interpolation
    if (linear_interpolation)
        if not(isempty(li_taps))
            if not(any(li_taps==1))
                % first reference
                H_CHE(:,1)=H_CHE(:,li_taps(1));
                li_taps=[1,li_taps];
            end
            if not(any(li_taps==pattern_len))
                % last reference
                H_CHE(:,pattern_len)=H_CHE(:,li_taps(end));
                li_taps=[li_taps,pattern_len];
            end
            for isym=1:(length(li_taps)-1)
                isym1=li_taps(isym);
                isym2=li_taps(isym+1);
                for isym_corr=(isym1+1):(isym2-1)
                    % H_est(k,i1p,:)+(H_est(k,i1m,:)-H_est(k,i1p,:))*(i1corr-i1m)/(i1p-i1m);
                    H_CHE(:,isym_corr)=H_CHE(:,isym1)+(H_CHE(:,isym2)-H_CHE(:,isym1))*(isym_corr-isym1)/(isym2-isym1);
                end
            end
        else
            warning('No pilots');
        end
    end
else
    warning('Incorrect size of SP.Pilots');
end
