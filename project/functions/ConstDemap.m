function [DATA,TYPE_DATA,TYPE_SYM] = ConstDemap(SYMB_stima,SP,idx_frame,idx_symb)

%-------------------------------------------------------------------------%
%                              ConstDemap.m
%-------------------------------------------------------------------------%
%   Function which demaps the symbols of the costellation ( QPSK, 16QAM or 
%   64QAM ) with Gray's encoding in bits associated to them. 
%   Each demodulator is defined in a subfunction.
% 
%   INPUT:
%    - SYMB_stima     Vector of symbols exstimated by the decisor
%    - SP             Simulator parameters (type: struct)
%    - idx_frame
%    - idx_symb
% 
%   OUTPUT:
%    - DATA           Vector of bits demapped from symbols
%    - TYPE_DATA, TYPE_SYM      1 - data, 2 - pilot
%
%   Autore: Puccitelli Marco
%-------------------------------------------------------------------------%

if (nargin==1)
    idx_frame=[];
    idx_symb=[];
end
if (nargin==2)
    idx_symb=[];
end

if not(isempty(idx_frame))&&not(isempty(idx_symb))
    % Check for pilots
    pattern_len=size(SP.Pilots,2);
    pilots_corr=SP.Pilots(:,rem((idx_frame-1)*SP.N_SCFDMAsymb+idx_symb-1,pattern_len)+1);
    i_pilots=find(pilots_corr');
end

% Inizialization of the vector of demapped bits
DATA = zeros(1,SP.Msymb*SP.mod_bitsymb);
TYPE_DATA = zeros(1,SP.Msymb*SP.mod_bitsymb);
TYPE_SYM = zeros(1,SP.Msymb);

for idx = 1:SP.Msymb 
    
    if not(any(i_pilots-idx))
        switch SP.Mod_type
            case 1
                data_tmp = QPSKdemod(SYMB_stima(idx));
            case 2
                data_tmp = QAM16demod(SYMB_stima(idx));
            case 3
                data_tmp = QAM64demod(SYMB_stima(idx));
            otherwise
                error('ConstDemap.m: Valore SP.Mod_type errato...');
        end
        type_data_tmp=ones(1,SP.mod_bitsymb);
    else
        data_tmp=NaN*ones(1,SP.mod_bitsymb);
        type_data_tmp=2*ones(1,SP.mod_bitsymb);
    end
    
    DATA((idx-1)*SP.mod_bitsymb+1:idx*SP.mod_bitsymb) = data_tmp;
    TYPE_DATA((idx-1)*SP.mod_bitsymb+1:idx*SP.mod_bitsymb) = type_data_tmp;
    TYPE_SYM(idx) = type_data_tmp(1);
    
end

%End_Of_Function
end


function bit = QPSKdemod(symb)
% Subfunction for QPSK demodulation
     
    switch symb 
        case 1/sqrt(2) + 1i* 1/sqrt(2), bit = [0 0];
        case 1/sqrt(2) - 1i* 1/sqrt(2), bit = [0 1];
        case -1/sqrt(2) + 1i* 1/sqrt(2), bit = [1 0];
        case -1/sqrt(2) - 1i* 1/sqrt(2), bit = [1 1];
        otherwise
            error('QPSKdemod subfunction: input errato');
    end
    
end

function bit = QAM16demod(symb)
% Subfunction for 16QAM demodulation

    switch symb 
        case 1/sqrt(10) + 1i* 1/sqrt(10), bit = [0 0 0 0];
        case 1/sqrt(10) + 1i* 3/sqrt(10), bit = [0 0 0 1];
        case 3/sqrt(10) + 1i* 1/sqrt(10), bit = [0 0 1 0];
        case 3/sqrt(10) + 1i* 3/sqrt(10), bit = [0 0 1 1];
        case 1/sqrt(10) - 1i* 1/sqrt(10), bit = [0 1 0 0];
        case 1/sqrt(10) - 1i* 3/sqrt(10), bit = [0 1 0 1];
        case 3/sqrt(10) - 1i* 1/sqrt(10), bit = [0 1 1 0];
        case 3/sqrt(10) - 1i* 3/sqrt(10), bit = [0 1 1 1];
        case -1/sqrt(10) + 1i* 1/sqrt(10), bit = [1 0 0 0];
        case -1/sqrt(10) + 1i* 3/sqrt(10), bit = [1 0 0 1];
        case -3/sqrt(10) + 1i* 1/sqrt(10), bit = [1 0 1 0];
        case -3/sqrt(10) + 1i* 3/sqrt(10), bit = [1 0 1 1];
        case -1/sqrt(10) - 1i* 1/sqrt(10), bit = [1 1 0 0];
        case -1/sqrt(10) - 1i* 3/sqrt(10), bit = [1 1 0 1];
        case -3/sqrt(10) - 1i* 1/sqrt(10), bit = [1 1 1 0];
        case -3/sqrt(10) - 1i* 3/sqrt(10), bit = [1 1 1 1];
        otherwise
            error('QAM16mod subfunction: input errato');
    end

end
 
function bit = QAM64demod(symb)
% Subfunction for 64QAM demodulation

    switch symb 
        case 3/sqrt(42) + 1i* 3/sqrt(42), bit = [0 0 0 0 0 0];
        case 3/sqrt(42) + 1i* 1/sqrt(42), bit = [0 0 0 0 0 1];
        case 1/sqrt(42) + 1i* 3/sqrt(42), bit = [0 0 0 0 1 0];
        case 1/sqrt(42) + 1i* 1/sqrt(42), bit = [0 0 0 0 1 1];
        case 3/sqrt(42) + 1i* 5/sqrt(42), bit = [0 0 0 1 0 0];  
        case 3/sqrt(42) + 1i* 7/sqrt(42), bit = [0 0 0 1 0 1];
        case 1/sqrt(42) + 1i* 5/sqrt(42), bit = [0 0 0 1 1 0];
        case 1/sqrt(42) + 1i* 7/sqrt(42), bit = [0 0 0 1 1 1];
        case 5/sqrt(42) + 1i* 3/sqrt(42), bit = [0 0 1 0 0 0];
        case 5/sqrt(42) + 1i* 1/sqrt(42), bit = [0 0 1 0 0 1];
        case 7/sqrt(42) + 1i* 3/sqrt(42), bit = [0 0 1 0 1 0];
        case 7/sqrt(42) + 1i* 1/sqrt(42), bit = [0 0 1 0 1 1];
        case 5/sqrt(42) + 1i* 5/sqrt(42), bit = [0 0 1 1 0 0];
        case 5/sqrt(42) + 1i* 7/sqrt(42), bit = [0 0 1 1 0 1];
        case 7/sqrt(42) + 1i* 5/sqrt(42), bit = [0 0 1 1 1 0];
        case 7/sqrt(42) + 1i* 7/sqrt(42), bit = [0 0 1 1 1 1];
        case 3/sqrt(42) - 1i* 3/sqrt(42), bit = [0 1 0 0 0 0];
        case 3/sqrt(42) - 1i* 1/sqrt(42), bit = [0 1 0 0 0 1];
        case 1/sqrt(42) - 1i* 3/sqrt(42), bit = [0 1 0 0 1 0];
        case 1/sqrt(42) - 1i* 1/sqrt(42), bit = [0 1 0 0 1 1];
        case 3/sqrt(42) - 1i* 5/sqrt(42), bit = [0 1 0 1 0 0];
        case 3/sqrt(42) - 1i* 7/sqrt(42), bit = [0 1 0 1 0 1];
        case 1/sqrt(42) - 1i* 5/sqrt(42), bit = [0 1 0 1 1 0];
        case 1/sqrt(42) - 1i* 7/sqrt(42), bit = [0 1 0 1 1 1];
        case 5/sqrt(42) - 1i* 3/sqrt(42), bit = [0 1 1 0 0 0];
        case 5/sqrt(42) - 1i* 1/sqrt(42), bit = [0 1 1 0 0 1];
        case 7/sqrt(42) - 1i* 3/sqrt(42), bit = [0 1 1 0 1 0];
        case 7/sqrt(42) - 1i* 1/sqrt(42), bit = [0 1 1 0 1 1];
        case 5/sqrt(42) - 1i* 5/sqrt(42), bit = [0 1 1 1 0 0];
        case 5/sqrt(42) - 1i* 7/sqrt(42), bit = [0 1 1 1 0 1];
        case 7/sqrt(42) - 1i* 5/sqrt(42), bit = [0 1 1 1 1 0];
        case 7/sqrt(42) - 1i* 7/sqrt(42), bit = [0 1 1 1 1 1];
        case -3/sqrt(42) + 1i* 3/sqrt(42), bit = [1 0 0 0 0 0];
        case -3/sqrt(42) + 1i* 1/sqrt(42), bit = [1 0 0 0 0 1];
        case -1/sqrt(42) + 1i* 3/sqrt(42), bit = [1 0 0 0 1 0];
        case -1/sqrt(42) + 1i* 1/sqrt(42), bit = [1 0 0 0 1 1];
        case -3/sqrt(42) + 1i* 5/sqrt(42), bit = [1 0 0 1 0 0];
        case -3/sqrt(42) + 1i* 7/sqrt(42), bit = [1 0 0 1 0 1];
        case -1/sqrt(42) + 1i* 5/sqrt(42), bit = [1 0 0 1 1 0];
        case -1/sqrt(42) + 1i* 7/sqrt(42), bit = [1 0 0 1 1 1];
        case -5/sqrt(42) + 1i* 3/sqrt(42), bit = [1 0 1 0 0 0];
        case -5/sqrt(42) + 1i* 1/sqrt(42), bit = [1 0 1 0 0 1];
        case -7/sqrt(42) + 1i* 3/sqrt(42), bit = [1 0 1 0 1 0];
        case -7/sqrt(42) + 1i* 1/sqrt(42), bit = [1 0 1 0 1 1];
        case -5/sqrt(42) + 1i* 5/sqrt(42), bit = [1 0 1 1 0 0];
        case -5/sqrt(42) + 1i* 7/sqrt(42), bit = [1 0 1 1 0 1];
        case -7/sqrt(42) + 1i* 5/sqrt(42), bit = [1 0 1 1 1 0];
        case -7/sqrt(42) + 1i* 7/sqrt(42), bit = [1 0 1 1 1 1];
        case -3/sqrt(42) - 1i* 3/sqrt(42), bit = [1 1 0 0 0 0];
        case -3/sqrt(42) - 1i* 1/sqrt(42), bit = [1 1 0 0 0 1];
        case -1/sqrt(42) - 1i* 3/sqrt(42), bit = [1 1 0 0 1 0];
        case -1/sqrt(42) - 1i* 1/sqrt(42), bit = [1 1 0 0 1 1];
        case -3/sqrt(42) - 1i* 5/sqrt(42), bit = [1 1 0 1 0 0];
        case -3/sqrt(42) - 1i* 7/sqrt(42), bit = [1 1 0 1 0 1];
        case -1/sqrt(42) - 1i* 5/sqrt(42), bit = [1 1 0 1 1 0];
        case -1/sqrt(42) - 1i* 7/sqrt(42), bit = [1 1 0 1 1 1];
        case -5/sqrt(42) - 1i* 3/sqrt(42), bit = [1 1 1 0 0 0];
        case -5/sqrt(42) - 1i* 1/sqrt(42), bit = [1 1 1 0 0 1];
        case -7/sqrt(42) - 1i* 3/sqrt(42), bit = [1 1 1 0 1 0];
        case -7/sqrt(42) - 1i* 1/sqrt(42), bit = [1 1 1 0 1 1];
        case -5/sqrt(42) - 1i* 5/sqrt(42), bit = [1 1 1 1 0 0];
        case -5/sqrt(42) - 1i* 7/sqrt(42), bit = [1 1 1 1 0 1];
        case -7/sqrt(42) - 1i* 5/sqrt(42), bit = [1 1 1 1 1 0];
        case -7/sqrt(42) - 1i* 7/sqrt(42), bit = [1 1 1 1 1 1]; 
        otherwise
            error('QAM64mod subfunction: input errato');
    end
    
end

%End_Of_File