function [tau,paths_OR,tau_sampling]=ext_channel(type,vel,f0,f_m,durata,OR)
%
% [tau,paths_OR,tau_sampling]=ext_channel(type,vel,f0,f_m,durata,OR)
%
% Extended channel models 3GPP 36.104
% EPA, EVA, ETU
% A multipath fading propagation condition is defined by a combination of a multi-path delay profile 
% and a maximum Doppler frequency fD which is either 5, 70 or 300 Hz. In addidion, 200 Hz Doppler frequency 
% is specified for UL timing adjustment performance requirement.
%
% INPUT
% type   = 'EPA', 'EVA, 'ETU'
% vel    = receiver speed while he's moving [km/h]
% f0     = carrier frequency [GHz]
% f_m    = maximum Doppler spread [if not NaN/[] -> "vel" not considerated]
% durata = total simulation duration [seconds]
% OR     = observation frequency [Hz], with resolution of 0.1,
%          resolution can be increased modifying "Dop_res" 
% OUTPUT
% tau    = vector of latency [nanoseconds]
% paths_OR = matrix that gives, for each latency, the sequence of complex
%            values of the response to the impulse.
%            size(path_OR)=[length(tau), durata*OR]

M=256;      % inputs' number of the Doppler filter (FIR)
Dop_res=1;  % Doppler resolution(Hz)
res_accu=20; % accuracy of resampling

%
switch type
    case 'EPA'
        PDP=[0.0,-1.0,-2.0,-3.0,-8.0,-17.2,-20.8];
        tau=[0,30,70,90,110,190,410];
        
    case 'EVA'
        PDP=[0.0,-1.5,-1.4,-3.6,-0.6,-9.1,-7.0,-12.0,-16.9];
        tau=[0,30,150,310,370,710,1090,1730,2510];
        
    case 'ETU'
        PDP=[-1.0,-1.0,-1.0,0.0,0.0,0.0,-3.0,-5.0,-7.0];
        tau=[0,50,120,200,230,500,1600,2300,5000];
        
    case 'RAYL'
        PDP=[0.0];
        tau=[0];

    otherwise
        error('Channel type unknown');
end

L=length(PDP);  %L number of taps

% calcolo massimo spostamento doppler f_m=(v/lambda)=(vel/3.6)*f0/c
c=3*10^8;
if (isnan(f_m))||(isempty(f_m))
    f_m=floor(vel*f0*10^9/(3.6*c*Dop_res))*Dop_res;    
end

% N = number of realizations with sampling rate SR = 2 * f_m
%SR=10
SR=floor(f_m*2/Dop_res)*Dop_res;
%
% N=floor(durata*SR); % to avoid errors
% N=max(1,floor(durata*SR));

if (OR<SR)
    % undersampling
    over_sampling=ceil((SR/OR));
    OR_eff=OR*over_sampling;
else
    over_sampling=1;
    OR_eff=OR;
end

%%% for the oversampliong (OR>SR) of the N realizations
if (SR>0)
    N=max(2,ceil(durata*SR));
    durata_new=N/SR;
    tau_sampling=1/OR;
    Nobs=max(1,floor(durata*OR));
    %
    m=lcm(round(SR/Dop_res),round(OR_eff/Dop_res));
    P=m/SR*Dop_res;
    Q=m/OR_eff*Dop_res;
    % P/Q = OR/SR
else
    % vel = 0
    N=1;
    durata_new=durata;
    tau_sampling=1/OR;
    Nobs=max(1,floor(durata*OR));
    %
    Q=1;
    P=1;
end

%%%%%%%%% generation of N independent realizations 
% (with Rayleigh distribution) 

P_=10.^(PDP/10); 
Fnorm=sum(P_); % normalization factor of the power
P_=P_/Fnorm;
%returns an L-by-N array of random integers
aux1=(randn(L,N)+1i*randn(L,N));
paths=sqrt(1/2)*aux1.*((sqrt(P_))'*ones(1,N));

% correlation between the N cahnnel's realizations with "classic" Doppler
% spectrum
if (M ~= 0)&&(N>1)
    for p=1:L
        D=1/2;    
        f_=[0:M*D]/(M*D);  % vector of normalized frequencies, between F=0 and F=1
        f1=f_(1:length(f_)-1);% vettore to avoid warning PSD(f_(end))=Inf !!!
        PSD=[1./sqrt(1-f1.^2),1./sqrt(1-f1(end)^2)];
        %%% instead of PSD=1./sqrt(1-f_.^2);   PSD(end)=PSD(end-1); 

        filt1=[PSD(1:end-1) zeros(1,M-floor(2*M*D)) PSD(end:-1:2)]; % spectrum S(f)
        % ripiegato e imbottito
        filt2=sqrt(filt1);% amplitude response of the filter 
        aux2=ifft(filt2);
        filt=ifftshift(aux2); % impulse response
        filt=real(filt); % elimination of "imaginary part" residuals
        filt=filt/sqrt(sum(filt.^2)); % normalization of the filter's energy
        % filtering the random sequence using the FIR with weigth 'filt'
        aux3=fftfilt(filt, [paths(p,:) zeros(1,M)]);
        paths_(p,:)=aux3(1+M/2:end-M/2);
    end
else
    paths_=paths;
end

% Upsampling Sequence  
if (N>1)
    for p=1:L
        paths_OR(p,:)=resample(paths_(p,:),round(P),round(Q),res_accu);
    end
    paths_OR=paths_OR(:,1:over_sampling:end);
    paths_OR=paths_OR(:,1:Nobs);
else
    paths_OR=paths_;
    paths_OR=repmat(paths_OR,1,Nobs);
end


