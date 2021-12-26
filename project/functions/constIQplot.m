function constIQplot(SP,Output)

%-------------------------------------------------------------------------%
%                             constIQplot.m
%-------------------------------------------------------------------------%
%  Function which represents the decision regions and the "scatterplot" IQ
%  of the received symbols. In addition, wrong symbols are reported.
% 
%  INPUT:
%   - SP             Simulator parameters (type: struct)
%   - Output         Vector of demodulated symbols (type: struct)  
% 
%  OUTPUT:
%   - "Scatter plot" of received symbols: wrong symbols are
%      highlighted by a black circle.
%     
%   Autore: Puccitelli Marco
%-------------------------------------------------------------------------%

switch SP.Mod_type
    case 1 
        % QPSK costellation
        costellazione = [ 1/sqrt(2) + 1i* 1/sqrt(2)
                          1/sqrt(2) - 1i* 1/sqrt(2)
                         -1/sqrt(2) + 1i* 1/sqrt(2)
                         -1/sqrt(2) - 1i* 1/sqrt(2) ];
                     
         plot(costellazione,'.'), hold on,
         line([-3 3],[0 0]), line([0 0],[-3 3]),
         plot(Output.SY_RX,'r.'),
         plot(Output.SY_RX(find(Output.SY_err)),'ko'),
         axis square, xlabel('I'), ylabel('Q'), 
         xlim([-2 2]), ylim([-2 2]),
         title(['Modulazione: ' SP.modulazione]);
                     
    case 2
        % 16QAM costellation
        costellazione = [ 1/sqrt(10) + 1i* 1/sqrt(10)
                          1/sqrt(10) + 1i* 3/sqrt(10)
                          3/sqrt(10) + 1i* 1/sqrt(10)
                          3/sqrt(10) + 1i* 3/sqrt(10)
                          1/sqrt(10) - 1i* 1/sqrt(10)
                          1/sqrt(10) - 1i* 3/sqrt(10)
                          3/sqrt(10) - 1i* 1/sqrt(10)
                          3/sqrt(10) - 1i* 3/sqrt(10)
                         -1/sqrt(10) + 1i* 1/sqrt(10)
                         -1/sqrt(10) + 1i* 3/sqrt(10)
                         -3/sqrt(10) + 1i* 1/sqrt(10)
                         -3/sqrt(10) + 1i* 3/sqrt(10)
                         -1/sqrt(10) - 1i* 1/sqrt(10)
                         -1/sqrt(10) - 1i* 3/sqrt(10)
                         -3/sqrt(10) - 1i* 1/sqrt(10)
                         -3/sqrt(10) - 1i* 3/sqrt(10) ];
                     
         voronoi(real(costellazione),imag(costellazione)), hold on,
         plot(Output.SY_RX,'r.'),
         plot(Output.SY_RX(find(Output.SY_err)),'ko'),
         axis square, xlabel('I'), ylabel('Q'), 
         xlim([-2 2]), ylim([-2 2]),
         title(['Modulazione: ' SP.modulazione]);

    case 3
        % 64QAM costellation
        costellazione = [ 3/sqrt(42) + 1i* 3/sqrt(42)
                          3/sqrt(42) + 1i* 1/sqrt(42)
                          1/sqrt(42) + 1i* 3/sqrt(42)
                          1/sqrt(42) + 1i* 1/sqrt(42)
                          3/sqrt(42) + 1i* 5/sqrt(42)  
                          3/sqrt(42) + 1i* 7/sqrt(42)
                          1/sqrt(42) + 1i* 5/sqrt(42)
                          1/sqrt(42) + 1i* 7/sqrt(42)
                          5/sqrt(42) + 1i* 3/sqrt(42)
                          5/sqrt(42) + 1i* 1/sqrt(42)
                          7/sqrt(42) + 1i* 3/sqrt(42)
                          7/sqrt(42) + 1i* 1/sqrt(42)
                          5/sqrt(42) + 1i* 5/sqrt(42)
                          5/sqrt(42) + 1i* 7/sqrt(42)
                          7/sqrt(42) + 1i* 5/sqrt(42)
                          7/sqrt(42) + 1i* 7/sqrt(42)
                          3/sqrt(42) - 1i* 3/sqrt(42)
                          3/sqrt(42) - 1i* 1/sqrt(42)
                          1/sqrt(42) - 1i* 3/sqrt(42)
                          1/sqrt(42) - 1i* 1/sqrt(42)
                          3/sqrt(42) - 1i* 5/sqrt(42)
                          3/sqrt(42) - 1i* 7/sqrt(42)
                          1/sqrt(42) - 1i* 5/sqrt(42)
                          1/sqrt(42) - 1i* 7/sqrt(42)
                          5/sqrt(42) - 1i* 3/sqrt(42)
                          5/sqrt(42) - 1i* 1/sqrt(42)
                          7/sqrt(42) - 1i* 3/sqrt(42)
                          7/sqrt(42) - 1i* 1/sqrt(42)
                          5/sqrt(42) - 1i* 5/sqrt(42)
                          5/sqrt(42) - 1i* 7/sqrt(42)
                          7/sqrt(42) - 1i* 5/sqrt(42)
                          7/sqrt(42) - 1i* 7/sqrt(42)
                         -3/sqrt(42) + 1i* 3/sqrt(42)
                         -3/sqrt(42) + 1i* 1/sqrt(42)
                         -1/sqrt(42) + 1i* 3/sqrt(42)
                         -1/sqrt(42) + 1i* 1/sqrt(42)
                         -3/sqrt(42) + 1i* 5/sqrt(42)
                         -3/sqrt(42) + 1i* 7/sqrt(42)
                         -1/sqrt(42) + 1i* 5/sqrt(42)
                         -1/sqrt(42) + 1i* 7/sqrt(42)
                         -5/sqrt(42) + 1i* 3/sqrt(42)
                         -5/sqrt(42) + 1i* 1/sqrt(42)
                         -7/sqrt(42) + 1i* 3/sqrt(42)
                         -7/sqrt(42) + 1i* 1/sqrt(42)
                         -5/sqrt(42) + 1i* 5/sqrt(42)
                         -5/sqrt(42) + 1i* 7/sqrt(42)
                         -7/sqrt(42) + 1i* 5/sqrt(42)
                         -7/sqrt(42) + 1i* 7/sqrt(42)
                         -3/sqrt(42) - 1i* 3/sqrt(42)
                         -3/sqrt(42) - 1i* 1/sqrt(42)
                         -1/sqrt(42) - 1i* 3/sqrt(42)
                         -1/sqrt(42) - 1i* 1/sqrt(42)
                         -3/sqrt(42) - 1i* 5/sqrt(42)
                         -3/sqrt(42) - 1i* 7/sqrt(42)
                         -1/sqrt(42) - 1i* 5/sqrt(42)
                         -1/sqrt(42) - 1i* 7/sqrt(42)
                         -5/sqrt(42) - 1i* 3/sqrt(42)
                         -5/sqrt(42) - 1i* 1/sqrt(42)
                         -7/sqrt(42) - 1i* 3/sqrt(42)
                         -7/sqrt(42) - 1i* 1/sqrt(42)
                         -5/sqrt(42) - 1i* 5/sqrt(42)
                         -5/sqrt(42) - 1i* 7/sqrt(42)
                         -7/sqrt(42) - 1i* 5/sqrt(42)
                         -7/sqrt(42) - 1i* 7/sqrt(42)];
                     
        voronoi(real(costellazione),imag(costellazione)), hold on,
        plot(Output.SY_RX,'r.'),
        plot(Output.SY_RX(find(Output.SY_err)),'ko'), 
        axis square, xlabel('I'), ylabel('Q'), 
        xlim([-2 2]), ylim([-2 2]),
        title(['Modulazione: ' SP.modulazione]);
                     
    otherwise
        error('Costellazione non definita!');      
        
end

%End_Of_Function
end

%End_Of_File