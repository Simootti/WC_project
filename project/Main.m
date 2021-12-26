clc
clear
close all

% Number of users N for each BS (that could also not be assigned to any BS)
N=10;
% Cell's radius
R=3e3;
% Vector representing the number of antennas for each Base Station
n_ant_RX_vect=[1 2 4 8];

addpath('./functions')
Inizializzazione_celle;

Sensitivity=-100; %dBm
shadowing = 'uniforme'; % or set as 'non_uniforme'
std_db=6; % standard deviation of the shadowing
f_c=5e9; % frequency carrier

Inizializzazione_Utenti;
Debug_interference;

%% Organize the User transmissions and the choice of the combiner
%% Run: Run_LTE_uplink


%% To plot the graph of the BER
%% Run: BER_MRC_MAX_SEL
