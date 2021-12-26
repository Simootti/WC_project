function [y,cp,len_cp] = RemCP(symbSCFDMA,idx_symb,SP) 
%
% [y,cp,len_cp] = RemCP(symbSCFDMA,idx_symb,SP)
%
%   Function that removes the cylic prefix from the received SCFDMA symbol.
%   In case of "Normal CP" (SP.CP_type=1), the length of the cyclic prefix
%   of the fisrt symbol of each slot (symbols 1 and 8 of the subframe) is 
%   bigger than the other symbols.  
% 
%   INPUT:
%    - symbSCFDMA     Vector of samples of received SCFDMA symbol
%    - idx_symb       Index of SCFDMA symbol
%    - SP             Simulator parameters (type: struct)
% 
%   OUTPUT:
%    - y              Vector of samples to be passed to DFT
%     
%   Autore: Puccitelli Marco
%-------------------------------------------------------------------------%

y=[];
cp=[];
len_cp=NaN;

switch SP.CP_type    
    
    case 1  % Normal CP 
        if (idx_symb == 1)
            len_cp=round((SP.CP(1)*SP.fs));
        elseif (idx_symb == 8)
            len_cp=round((SP.CP(1)*SP.fs));
        else            
            len_cp=round((SP.CP(2)*SP.fs));
        end  
        y = symbSCFDMA(len_cp+1:end);
        cp = symbSCFDMA(1:len_cp);
        
    case 2  % Extended CP
        len_cp=round((SP.CP(1)*SP.fs));
        y = symbSCFDMA(len_cp+1:end);
        cp = symbSCFDMA(1:len_cp);
              
    otherwise
        error('RemCP.m: valore SP.CP_type errato...');
end

%End_Of_Function
end

%End_Of_File