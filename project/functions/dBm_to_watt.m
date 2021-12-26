%% Function to convert dBm to Watt

function [watt_power] = dBm_to_watt(dBm_power)

    % returns the linear version of the dBm power
    watt_power = (10^(dBm_power/10))/1000;
    
end