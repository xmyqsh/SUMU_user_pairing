clear all
clc

global user_selected channel

s=rand('state');
sn=randn('state');
ctrlpar=ctrlparset();
ctrlpar.ntones = 1;
ctrlpar.ntimesamples = 1;
ctrlpar.FBgran = 'SB Level';    %'RB Level';

ctrlpar.EESMbetas = [5,5.01,5.01,0.84,1.67,1.61,1.64,3.87,5.06,6.4,12.59,17.59,23.33,29.45,33.05,35.41];  % weigthed EESM beta values
ctrlpar.MCSs = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
ctrlpar.TargetSINR = [-7.5176,-5.6776,-3.5476,-1.5076,0.4924,2.3924,4.3972,6.2272,8.1972,10.1191,12.0191,13.9691,15.8191,17.7491,19.6391];


ctrlpar.nusers = 25;
ctrlpar.ndrops = 1000;  %100;
ctrlpar.Nt = 8;
ctrlpar.Nr = 2;
%   antenna error
ctrlpar.phase_error = [-5, 5];  %   degree
ctrlpar.amp_error = [-1, 1];    %   dB

%----------------------%
% capi=sim_MUCB(ctrlpar);      % 5bit Feedback No AutoSwitch

eval(['load channel_', num2str(ctrlpar.nusers), 'u_8t_2r.mat'])
eval(['load user_selected_', num2str(ctrlpar.nusers-20), '.mat'])

eval(['user_selected = user_selected_', num2str(ctrlpar.nusers-20), ';']);
ebno_range = ctrlpar.ebno_range;

for ebno_ind = 1:length(ebno_range)
    ctrlpar.ebno_range = ebno_range(ebno_ind);
    capi=sim_MUCB_selected_user(ctrlpar);
end
