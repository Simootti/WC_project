
function [SYMB,reference_symbol,pilots_corr] = ConstMap(DATA,SP,idx_frame,idx_symb)

%-------------------------------------------------------------------------%
%                              ConstMap.m
%-------------------------------------------------------------------------%
%   Function which maps the bits that must be transmitted  
%   in the symbols of the chosen costellation
%   ( QPSK, 16QAM or 64QAM ) applying the Gray encoding. 
%   Each modulator is defined in a subfunction.
%   Checks and adds the "reference symbols / pilots"
% 
%   INPUT:
%    - DATA           Vector of bits to be mapped in the costellation
%    - SP             Simulator parameters (type: struct)
%    - idx_frame
%    - idx_symb
%
%   OUTPUT:
%    - SYMB           Vector of complex symbols
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

% Disposition of bits related to symbols in columns.
% "reshape" makes DATA becaome a matrix 
% of "SP.mod_bitsymb * SP.Msymb" elements
% SP.mod_bitsymb
% SP.Msymb
tmp = reshape(DATA,SP.mod_bitsymb,SP.Msymb);
reference_symbol=0;
pilots_corr=[];

for idx = 1:SP.Msymb     
    switch SP.Mod_type
        case 1
            SYMB(idx) = QPSKmod(tmp(:,idx)'); 
        case 2
            SYMB(idx) = QAM16mod(tmp(:,idx)'); 
        case 3
            SYMB(idx) = QAM64mod(tmp(:,idx)');
        otherwise
            error('ConstMap.m: Valore SP.Mod_type errato...');
    end  
end

if not(isempty(idx_frame))&&not(isempty(idx_symb))
    % Check for pilots
    pattern_len=size(SP.Pilots,2);
    ip_corr=rem((idx_frame-1)*SP.N_SCFDMAsymb+idx_symb-1,pattern_len)+1;
    flag_pilots_corr=SP.Flag_Pilots(:,ip_corr);
    i_pilots=find(flag_pilots_corr.');
    pilots_corr=SP.Pilots(i_pilots,ip_corr);
    if not(isempty(i_pilots))
        SYMB(i_pilots)=pilots_corr.';
        if (length(i_pilots)==SP.Msymb)
            reference_symbol=1;
        end 
    end
end

%End_Of_Function
end


function symb = QPSKmod(bit)
% Subfunction for QPSK modulation
 
    s = strcat(num2str(bit(1)),num2str(bit(2)));
    switch s 
        case '00', symb = 1/sqrt(2) + 1i* 1/sqrt(2);
        case '01', symb = 1/sqrt(2) - 1i* 1/sqrt(2);
        case '10', symb = -1/sqrt(2) + 1i* 1/sqrt(2);
        case '11', symb = -1/sqrt(2) - 1i* 1/sqrt(2);
        otherwise
            error('QPSKmod subfunction: input errato');
    end
    
end

function symb = QAM16mod(bit)
% Subfunction for 16QAM modulation

    s = strcat(num2str(bit(1)),num2str(bit(2)), ... 
               num2str(bit(3)),num2str(bit(4)));
    switch s 
        case '0000', symb = 1/sqrt(10) + 1i* 1/sqrt(10);
        case '0001', symb = 1/sqrt(10) + 1i* 3/sqrt(10);
        case '0010', symb = 3/sqrt(10) + 1i* 1/sqrt(10);
        case '0011', symb = 3/sqrt(10) + 1i* 3/sqrt(10);
        case '0100', symb = 1/sqrt(10) - 1i* 1/sqrt(10);
        case '0101', symb = 1/sqrt(10) - 1i* 3/sqrt(10);
        case '0110', symb = 3/sqrt(10) - 1i* 1/sqrt(10);
        case '0111', symb = 3/sqrt(10) - 1i* 3/sqrt(10);
        case '1000', symb = -1/sqrt(10) + 1i* 1/sqrt(10);
        case '1001', symb = -1/sqrt(10) + 1i* 3/sqrt(10);
        case '1010', symb = -3/sqrt(10) + 1i* 1/sqrt(10);
        case '1011', symb = -3/sqrt(10) + 1i* 3/sqrt(10);
        case '1100', symb = -1/sqrt(10) - 1i* 1/sqrt(10);
        case '1101', symb = -1/sqrt(10) - 1i* 3/sqrt(10);
        case '1110', symb = -3/sqrt(10) - 1i* 1/sqrt(10);
        case '1111', symb = -3/sqrt(10) - 1i* 3/sqrt(10);
        otherwise
            error('QAM16mod subfunction: input errato');
    end

end
 
function symb = QAM64mod(bit)
% Subfunction for 64QAM modulation

    s = strcat(num2str(bit(1)),num2str(bit(2)), ... 
               num2str(bit(3)),num2str(bit(4)), ...
               num2str(bit(5)),num2str(bit(6)));
    switch s 
        case '000000', symb = 3/sqrt(42) + 1i* 3/sqrt(42);
        case '000001', symb = 3/sqrt(42) + 1i* 1/sqrt(42);
        case '000010', symb = 1/sqrt(42) + 1i* 3/sqrt(42);
        case '000011', symb = 1/sqrt(42) + 1i* 1/sqrt(42);
        case '000100', symb = 3/sqrt(42) + 1i* 5/sqrt(42);  
        case '000101', symb = 3/sqrt(42) + 1i* 7/sqrt(42);
        case '000110', symb = 1/sqrt(42) + 1i* 5/sqrt(42);
        case '000111', symb = 1/sqrt(42) + 1i* 7/sqrt(42);
        case '001000', symb = 5/sqrt(42) + 1i* 3/sqrt(42);
        case '001001', symb = 5/sqrt(42) + 1i* 1/sqrt(42);
        case '001010', symb = 7/sqrt(42) + 1i* 3/sqrt(42);
        case '001011', symb = 7/sqrt(42) + 1i* 1/sqrt(42);
        case '001100', symb = 5/sqrt(42) + 1i* 5/sqrt(42);
        case '001101', symb = 5/sqrt(42) + 1i* 7/sqrt(42);
        case '001110', symb = 7/sqrt(42) + 1i* 5/sqrt(42);
        case '001111', symb = 7/sqrt(42) + 1i* 7/sqrt(42);
        case '010000', symb = 3/sqrt(42) - 1i* 3/sqrt(42);
        case '010001', symb = 3/sqrt(42) - 1i* 1/sqrt(42);
        case '010010', symb = 1/sqrt(42) - 1i* 3/sqrt(42);
        case '010011', symb = 1/sqrt(42) - 1i* 1/sqrt(42);
        case '010100', symb = 3/sqrt(42) - 1i* 5/sqrt(42);
        case '010101', symb = 3/sqrt(42) - 1i* 7/sqrt(42);
        case '010110', symb = 1/sqrt(42) - 1i* 5/sqrt(42);
        case '010111', symb = 1/sqrt(42) - 1i* 7/sqrt(42);
        case '011000', symb = 5/sqrt(42) - 1i* 3/sqrt(42);
        case '011001', symb = 5/sqrt(42) - 1i* 1/sqrt(42);
        case '011010', symb = 7/sqrt(42) - 1i* 3/sqrt(42);
        case '011011', symb = 7/sqrt(42) - 1i* 1/sqrt(42);
        case '011100', symb = 5/sqrt(42) - 1i* 5/sqrt(42);
        case '011101', symb = 5/sqrt(42) - 1i* 7/sqrt(42);
        case '011110', symb = 7/sqrt(42) - 1i* 5/sqrt(42);
        case '011111', symb = 7/sqrt(42) - 1i* 7/sqrt(42);
        case '100000', symb = -3/sqrt(42) + 1i* 3/sqrt(42);
        case '100001', symb = -3/sqrt(42) + 1i* 1/sqrt(42);
        case '100010', symb = -1/sqrt(42) + 1i* 3/sqrt(42);
        case '100011', symb = -1/sqrt(42) + 1i* 1/sqrt(42);
        case '100100', symb = -3/sqrt(42) + 1i* 5/sqrt(42);
        case '100101', symb = -3/sqrt(42) + 1i* 7/sqrt(42);
        case '100110', symb = -1/sqrt(42) + 1i* 5/sqrt(42);
        case '100111', symb = -1/sqrt(42) + 1i* 7/sqrt(42);
        case '101000', symb = -5/sqrt(42) + 1i* 3/sqrt(42);
        case '101001', symb = -5/sqrt(42) + 1i* 1/sqrt(42);
        case '101010', symb = -7/sqrt(42) + 1i* 3/sqrt(42);
        case '101011', symb = -7/sqrt(42) + 1i* 1/sqrt(42);
        case '101100', symb = -5/sqrt(42) + 1i* 5/sqrt(42);
        case '101101', symb = -5/sqrt(42) + 1i* 7/sqrt(42);
        case '101110', symb = -7/sqrt(42) + 1i* 5/sqrt(42);
        case '101111', symb = -7/sqrt(42) + 1i* 7/sqrt(42);
        case '110000', symb = -3/sqrt(42) - 1i* 3/sqrt(42);
        case '110001', symb = -3/sqrt(42) - 1i* 1/sqrt(42);
        case '110010', symb = -1/sqrt(42) - 1i* 3/sqrt(42);
        case '110011', symb = -1/sqrt(42) - 1i* 1/sqrt(42);
        case '110100', symb = -3/sqrt(42) - 1i* 5/sqrt(42);
        case '110101', symb = -3/sqrt(42) - 1i* 7/sqrt(42);
        case '110110', symb = -1/sqrt(42) - 1i* 5/sqrt(42);
        case '110111', symb = -1/sqrt(42) - 1i* 7/sqrt(42);
        case '111000', symb = -5/sqrt(42) - 1i* 3/sqrt(42);
        case '111001', symb = -5/sqrt(42) - 1i* 1/sqrt(42);
        case '111010', symb = -7/sqrt(42) - 1i* 3/sqrt(42);
        case '111011', symb = -7/sqrt(42) - 1i* 1/sqrt(42);
        case '111100', symb = -5/sqrt(42) - 1i* 5/sqrt(42);
        case '111101', symb = -5/sqrt(42) - 1i* 7/sqrt(42);
        case '111110', symb = -7/sqrt(42) - 1i* 5/sqrt(42);
        case '111111', symb = -7/sqrt(42) - 1i* 7/sqrt(42);
        otherwise
            error('QAM64mod subfunction: input errato');
    end
    
end

%End_Of_File