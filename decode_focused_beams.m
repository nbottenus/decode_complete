function [rf_fsa] = decode_focused_beams(rf,transmit_delays,method)
%DECODE_FOCUSED_BEAMS - Decode RF channel data from a set of focused 
% transmit beams into the complete data set using the applied focal timings
%
% rf_fsa = decode_focused_beams(rf,transmit_delays);
% rf_fsa = decode_focused_beams(rf,transmit_delays, METHOD);
%
% Decoding can be performed in either the time or frequency domains. The
% two methods produce the same result, but frequency is faster.
%
% Parameters:
% rf - RF data (time sample x receive channel x transmit event)
% transmit_delays - Transmit focal delays in samples (transmit event x transmit element)
% method - 'frequency' or 'time' (default: frequency)
%
% Returns:
% rf_fsa - Decoded complete data set (time sample x receive channel x transmit element)
%
% Author: Nick Bottenus
% Contact: nick.bottenus@duke.edu

% Check inputs
if(~exist('method','var'))
    method='frequency';
else
    assert(strcmp(method,'frequency')||strcmp(method,'time'),'Method must be ''frequency'' or ''time''')
end
[n_samples, n_receives, n_transmits]=size(rf);
n_elements=size(transmit_delays,2);
assert(size(transmit_delays,1)==n_transmits,'Transmit event dimensions inconsistent between rf and transmit_timings')

switch method
    
    % ======================================Frequency domain implementation
    case 'frequency'
       
    % 1-D FFT to convert time to frequency
    RF=fft(single(rf));
    RF=permute(RF,[2 3 1]); % (receive channel x transmit event x time sample)
    frequency=(0:n_samples-1)/n_samples;
    
    % Apply decoding matrix at each frequency
    RF_adj=zeros(n_samples,n_receives,n_elements,'single');
    for i=1:ceil(n_samples/2) % only compute half, assume symmetry
        omega=2*pi*frequency(i);
        Hinv=exp(1j*omega*transmit_delays);
        RF_adj(i,:,:)=RF(:,:,i)*Hinv;
    end
    
    % Inverse FFT for real signal
    rf_fsa=ifft(RF_adj,'symmetric');
    
    % ===========================================Time domain implementation
    case 'time'

    rf_fsa=zeros(n_samples,n_receives,n_elements,'single');
    samples=1:n_samples;
    % Interpolate the data for each transmit event with delays
    % corresponding to each recovered transmit element
    for j=1:n_transmits
        rf_fsa=rf_fsa+interp1(samples,single(rf(:,:,j)),...
            repmat(samples(:),1,n_elements)+repmat(transmit_delays(j,:),n_samples,1),...
            'linear');
    end
    rf_fsa(isnan(rf_fsa))=0;
end