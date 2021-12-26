function DATA = DataGen(SP) 

%-------------------------------------------------------------------------%
%                              DataGen.m
%-------------------------------------------------------------------------%
%   Function which allows to generate (SP.Data_type = 1) or read from file 
%   (SP.Data_type = 2) a sequence of bits.       
%   The returned number of bits is equal to the number of bits transmitted   
%   durinng the SCFDMA-symbol-time for the chosen modulation.
% 
%   INPUT:
%    - SP             Simulator parameters (type: struct)
% 
%   OUTPUT:
%    - DATA           Vector of data bits
%     
%   Autore: Puccitelli Marco
%-------------------------------------------------------------------------%

    switch SP.Data_type

        case 1
             rng('shuffle')
             DATA = round(rand(SP.Msymb*SP.mod_bitsymb,SP.N_SCFDMAsymb,'double'));

        case 2
            % Reading data from file
            error('Lettura bit da file non implementata ...')

        otherwise
            error('Valore Data_type errato!')
    end
    
end