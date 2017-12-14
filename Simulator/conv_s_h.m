function y=conv_s_h(s,h,pdp,Nc,samp_freq,tg)
%tg=5.21e-6;                       % guerad time -> 80 samples

delays=pdp(:,1);
delta_t=1/samp_freq;
Npaths = length(delays);                      % No. of paths considered for the channel
deltans=delta_t/1e-9;             % Sampling interval in ns


delays = round(delays/(deltans))+1;
Ng= round(tg/(delta_t))+1;
f_zeros=find(h==0);
h(f_zeros)=[];

aux = zeros(Npaths,Nc+Ng);

for n=1:Npaths
    

conv_sh=h(n)*s;

aux(n,delays(n):Nc+Ng-1+delays(n)-1)=conv_sh;

end


y=sum(aux);
y=y(1:Nc+Ng-1);