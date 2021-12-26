%% Pathloss function

function PL = PathLoss(d,fc,B)

   speedOfLight = 299792458.0; % Speed of light in m/s
   
   % Pathloss calculation
   PL=(4*pi*fc*d/speedOfLight)^2;
   
end