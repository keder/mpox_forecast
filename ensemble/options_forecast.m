 
function [getperformance, deletetempfiles, forecastingperiod, weight_type1]=options_forecast

% <==============================================================================>
% <========================== Forecasting parameters ===================================>
% <==============================================================================>

getperformance=1; % flag or indicator variable (1/0) to calculate forecasting performance metrics or not

deletetempfiles=1; %flag or indicator variable (1/0) to indicate whether we wan to delete Forecast..mat files after use

forecastingperiod=1; % forecast horizon (number of time units ahead)

% <==============================================================================>
% <====================== weighting scheme for ensemble model ============================>
% <==============================================================================>

weight_type1=1; % -1= equally weighted from the top models, 0= weighted ensemble based on AICc, 1= weighted ensemble based on relative likelihood (Akaike weights), 
% 2=weighted ensemble based on the weighted interval score of the calibration period (WISC).
