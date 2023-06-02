function SDR=eval_sdr2(Est_s,origin,fs,q)

% Est_s is the estimated signal (LxI or Lx1)
% origin is the target signal before mixing (LxI2 or Lx1)
% The output is the Signal-to-Distortion ratio

% Modified on 04/11/2010 15:17 by Atiyeh.A
% modified again by Qingju Liu
if nargin>2 && ~isempty(fs)
else
fs = 16000;
end
if nargin<4
q=0.032*fs; % a Wiener filter of 32 ms
end

if isvector(Est_s)
Est_s = Est_s(:);
else
    L1 = size(Est_s,1); L2 = size(Est_s,2); 
    if L1<L2
        Est_s = Est_s.';
    end
end
[L,I] = size(Est_s);

if isvector(origin)
origin = origin(:);
else
    L1 = size(origin,1); L2 = size(origin,2); 
    if L1<L2
        origin = origin.';
    end
end
s=origin(1:L,:);% Because reconstruct function changes the length

for i = 1:I
    Estim_s_temp=Est_s(:,i);
    % Any energy in the estimated signal that can be explained with a
    % linear combination of delayed versions of the target signal
    s1 = s(:,min(i,size(s,2)));
    h_wien=wiener_filter2(s1,Estim_s_temp,q);
    S1_hat=filter(h_wien,1,s1);
    
    target_energy=var(S1_hat);
    distortion_energy=var(Estim_s_temp-S1_hat);
    
    SDR(i)=10*log10(target_energy/distortion_energy);
end

SDR = mean(SDR);
