%-------------------------------------------------------------------------%
%                              AddCP.m
%-------------------------------------------------------------------------%
%   Function which adds the cyclic prefix for the transmission of the
%   SCFDMA symbols. In case of Normal CP (SP.CP_type=1), the length of
%   the cyclic prefix of each slot's first symbol (symbols 1 and 8 of the  
%   subframe) is bigger than the other symbols.  
% 
%   INPUT:  
%    - y              "Samples vector" of the IDFT
%    - idx_symb       Index of the SCFDMA symbol
%    - SP             Simulator's parameters (type: struct)
% 
%   OUTPUT:
%    - symbSCFDMA      "Samples vector" to transmit
%     
%   Autore: Puccitelli Marco
%-------------------------------------------------------------------------%

function symbSCFDMA = AddCP(y,idx_symb,SP) 


switch SP.CP_type    
    
    case 1  % Normal Cyclic Prefix 
        if (idx_symb == 1)
            % If the symbol is the number 1-st (or 8-th), we insert the
            % cyclic prefix because in the Normal Cyclic Prefix each  
            % Time Slot has 7 OFDMA symbols (only 6 in the Extended CP case)
            % From 
            % ""SimConfig.m:{
            % SP.CP(1)=160*SP.Ts, Ts:basic time unit, lunghezza cyclic
            % prefix
            % SP.fs=frequenza di campionamento
            % }""
            
            % "round" approximates to the nearest integer
            % "y(end - round(..) : end)" takes the elements with indexes 
            % that will go from the difference "end - round(..)" till "end"
            % So "cp" is a vector with first index "end-round(..)+1" 
            % so that it will never be 0
            cp = y(end-round((SP.CP(1)*SP.fs))+1:end);
            % The SCFDMA symbol depends on the cyclic prefix  
            % and on the symbol resulting from the IFFT/IDFT 
            symbSCFDMA = [cp y];
        elseif (idx_symb == 8)        % same as before for the 8-th symbol
            cp = y(end-round((SP.CP(1)*SP.fs))+1:end);
            symbSCFDMA = [cp y];
        else
            % SP.CP(2)=144*SP.Ts < SP.CP(1)
            cp = y(end-round((SP.CP(2)*SP.fs))+1:end);
            symbSCFDMA = [cp y];
        end  
        
    case 2  % Extended CP
        % values of the cyclic prefix are the same for all the symbols
        cp = y(end-round((SP.CP(1)*SP.fs))+1:end);
        symbSCFDMA = [cp y];
              
    otherwise
        error('AddCp.m: valore SP.CP_type errato...');
end

%End_Of_Function
end

%End_Of_File