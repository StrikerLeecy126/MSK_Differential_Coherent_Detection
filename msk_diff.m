%#ok<*AGROW> 
% mskmod local data with mskmod function
Local_data=mskmod(Local_data_bit,nsamp,'nondiff');

% Add correlation function with different k
R1=zeros(1,length(I_data_LPF));

for k=1:2:bps
    [R,D,E]=diff_coh(Idata_LPF,Qdata_LPF,Local_data,Tb,fs,k);
    R1=R1+R;
end


% Tb: Real Time lasted in analogue for 1-bit.
% fs: Sampling Frequency

function [R,D,E]=diff_coh(Idata,Qdata,Local_data,Tb,fs,k)

% Local Baseband Data is a Complex Signal
Local_data_R=real(Local_data);
Local_data_I=imag(Local_data);

sps=Tb*fs;

% k-bit differential demodulation for received signal
D=[];
for m=1:length(Idata)-k*sps
    D(m)=Qdata(m+k*sps).*Idata(m)-Idata(m+k*sps).*Qdata(m);
end

% k-bit differential demodulation for local signal
E=[];
for m=1:length(Header_MSK)-k*sps
    E(m)=Local_data_I(m+k*sps).*Local_data_R(m)-Local_data_R(m+k*sps).*Local_data_I(m);
end

% Sliding Window Corrrlation
coh_window=length(E);
R=[];
for m=1:length(D)-coh_window
    R(m)=sum(D(m:m+coh_window-1).*E);
end

if(length(R)<length(Idata))
    R=[R,zeros(1,length(Idata)-length(R))];
end

end
