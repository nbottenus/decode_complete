% Load the sample data and transmit delay timings
load('sample_data.mat');

% Optionally reduce the number of transmit events in the saved data
% ds=6;
% rf=rf(:,:,1:ds:end);
% transmit_delays=transmit_delays(1:ds:end,:);

% Perform decoding to produce the complete data set
rf_fsa=decode_focused_beams(rf,transmit_delays);

% Perform a basic diverging wave focusing of the complete data set
x=linspace(-15,15,200)/1000;
z=linspace(15,60,500)/1000;
r=(params.t0+(0:size(rf_fsa,1)-1))/params.fs*params.c;
rf_focused=beamform(rf_fsa,r,params.rx_pos,x,z);

% Apply a high-pass filter to the data to remove interpolation artifacts
[b,a]=butter(2,500e3/(params.fs/2),'high');
rf_focused=filter(b,a,rf_focused);

% Display image
env=abs(hilbert(rf_focused));
env=env/max(env(:));
imagesc(x*1e3,z*1e3,db(env),[-50 0]);axis image;colormap gray
xlabel('Lateral (mm)')
ylabel('Axial (mm)')
title('Recovered complete data set')

% Load the saved figures
% openfig('dynamic_receive.fig')
% openfig('recovered_complete.fig')