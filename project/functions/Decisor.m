function SYMB_stima = Decisor(SYMB_RX,SP)

%-------------------------------------------------------------------------%
%                              Decisor.m
%-------------------------------------------------------------------------%
%   Function which returns the estimation of received symbols with an 
%   'Hard-Decision' method (symbol compared all the possible symbols:
%   then is chosen the one with smaller Hamming distance Hamming w.r.t
%   the considered symbol). {'Soft Decision'-> Euclidean Distance}
% 
%   INPUT:
%    - SYMB_RX        Vector of complex symbols received
%    - SP             Simulator parameters (type: struct)
% 
%   OUTPUT:
%    - SYMB_stima     Vector of symbols estimated by the decisor
%     
%   Autore: Puccitelli Marco
%-------------------------------------------------------------------------%

switch SP.Mod_type
    case 1 
        % QPSK decisor
        for idx_s = 1:length(SYMB_RX) 

            I = real(SYMB_RX(idx_s));
            Q = imag(SYMB_RX(idx_s));

            if (I >= 0)
                tmpI = 1/sqrt(2);       
            elseif (I < 0)
                tmpI = -1/sqrt(2);    
            end

            if (Q >= 0)
                tmpQ = 1/sqrt(2);                      
            elseif (Q < 0)
                tmpQ = -1/sqrt(2);    
            end
                
            SYMB_stima(idx_s) = tmpI + 1i*tmpQ;
        end
                     
    case 2
        % 16QAM decisor
        for idx_s = 1:length(SYMB_RX) 

            I = real(SYMB_RX(idx_s));
            Q = imag(SYMB_RX(idx_s));

            if (I >= 0 && I < 2/sqrt(10))
                tmpI = 1/sqrt(10);            
            elseif (I >= 2/sqrt(10))
                tmpI = 3/sqrt(10);
            elseif (I >= -2/sqrt(10) && I < 0)
                tmpI = -1/sqrt(10);            
            elseif (I < -2/sqrt(10))
                tmpI = -3/sqrt(10);    
            end

            if (Q >= 0 && Q < 2/sqrt(10))
                tmpQ = 1/sqrt(10);            
            elseif (Q >= 2/sqrt(10))
                tmpQ = 3/sqrt(10);
            elseif (Q >= -2/sqrt(10) && Q < 0)
                tmpQ = -1/sqrt(10);            
            elseif (Q < -2/sqrt(10))
                tmpQ = -3/sqrt(10);    
            end
                
            SYMB_stima(idx_s) = tmpI + 1i*tmpQ;
        end

    case 3   
        % 64QAM decisor
        for idx_s = 1:length(SYMB_RX) 

            I = real(SYMB_RX(idx_s));
            Q = imag(SYMB_RX(idx_s));

            if (I >= 0 && I < 2/sqrt(42))
                tmpI = 1/sqrt(42);
            elseif (I >= 2/sqrt(42) && I < 4/sqrt(42))
                tmpI = 3/sqrt(42);
            elseif (I >= 4/sqrt(42) && I < 6/sqrt(42))
                tmpI = 5/sqrt(42);
            elseif (I >= 6/sqrt(42))
                tmpI = 7/sqrt(42);
            elseif (I >= -2/sqrt(42) && I < 0)
                tmpI = -1/sqrt(42);
            elseif (I >= -4/sqrt(42) && I < -2/sqrt(42))
                tmpI = -3/sqrt(42);
            elseif (I >= -6/sqrt(42) && I < -4/sqrt(42))
                tmpI = -5/sqrt(42);
            elseif (I < -6/sqrt(42))
                tmpI = -7/sqrt(42);    
            end

            if (Q >= 0 && Q < 2/sqrt(42))
                tmpQ = 1/sqrt(42);
            elseif (Q >= 2/sqrt(42) && Q < 4/sqrt(42))
                tmpQ = 3/sqrt(42);
            elseif (Q >= 4/sqrt(42) && Q < 6/sqrt(42))
                tmpQ = 5/sqrt(42);
            elseif (Q >= 6/sqrt(42))
                tmpQ = 7/sqrt(42);
            elseif (Q >= -2/sqrt(42) && Q < 0)
                tmpQ = -1/sqrt(42);
            elseif (Q >= -4/sqrt(42) && Q < -2/sqrt(42))
                tmpQ = -3/sqrt(42);
            elseif (Q >= -6/sqrt(42) && Q < -4/sqrt(42))
                tmpQ = -5/sqrt(42);
            elseif (Q < -6/sqrt(42))
                tmpQ = -7/sqrt(42);    
            end
                
            SYMB_stima(idx_s) = tmpI + 1i*tmpQ;
        end
        
    otherwise
        error('Decisor.m: valore SP.Mod_type errato...');      
        
end

%End_Of_Function
end

%End_Of_File