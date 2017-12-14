%*****************Channel Generation *****************
function [ht, Hf]=channel_gen(pdp,samp_freq, Nc)

delta_t=1/samp_freq;              %sample duratiion
Npaths = length(pdp(:,1));        % No. of paths considered for the channel
deltans=delta_t/1e-9;             % Sampling interval in ns

path_pot_lin=10.^(pdp(:,2)/10);
path_pot_lin=path_pot_lin./sum(path_pot_lin);

delays = pdp(:,1);
delays = round(delays./(deltans))+1;

multipath = zeros(1,Npaths);

for n=1:Npaths
pot_componente=0.5;
multipath(n)=sqrt(pot_componente)*randn(1,1)+j*sqrt(pot_componente)*randn(1,1);
multipath(n)=multipath(n).*sqrt(path_pot_lin(n));

end

RI=zeros(1,Nc);
RI(delays) = RI(delays) + multipath;

ht=RI;
Hf=fft(ht);