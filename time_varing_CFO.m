% Carrier frequency offset recovery (PLL based) with system object and 
% a replica function. Much clear for reading, and the equations 
% can be found from Matlab docs.

% by Xiao Liu, Dalhousie University, UW-Stream Lab
% May 4th (be with you), 2018
%% create symbols with time varing phase offset

close all
clc
M = 4;
rng(1);
symbNum = 5e3;
msg = randi([0 M-1], symbNum, 1);
modData = pskmod(msg, M, pi/M);
t = (1:symbNum)';
phaseOff = pi*( 5e-3*t+sin(2e-3*t).*(1e-3+5e-4.*t));
% phaseOff = 0;

plot(t, phaseOff/pi,'.')
title('Phase offset')
noisyData = single(awgn(modData.*exp(1j*phaseOff), 10));
scatterplot(noisyData(1:end))
title('constellation with CFO')
% modDataIR = single([real(modData), imag(modData)]);
% sigRI = [real(noisyData), imag(noisyData)];
% sigPhase = angle(noisyData);
% save('symbCFO.mat','sigRI','sigPhase','modDataIR')
%% carrier offset recovery
% build-in system object fot carrier synchronization
carrierSync = comm.CarrierSynchronizer( ...
    'SamplesPerSymbol',1,...
    'Modulation','QPSK');
syncSignal = carrierSync(noisyData);
scatterplot(syncSignal(end-500:end))
title('carrier recovery with build-in system object')

[mySyncSig, err, phi, lambda] = carrSync(noisyData);
scatterplot(mySyncSig(end-500:end))
title('carrier recovery with a replica function')

figure
subplot(3,1,1)
plot(err)
title('error')
subplot(3,1,2)
plot(phi)
title('\phi')
subplot(3,1,3)
plot(lambda)
title('\lambda')
%% 
function [y, err, phi, lambda] = carrSync(x)
    % carrier sync for BPSK, QPSK
    % a replica of Matlab system object
    len = length(x);
    y = x;
    err = zeros(1,len);
    Bn = 0.01; % NormalizedLoopBandwidth
    xi = 0.707; % DampingFactor
    Kp = 2; % for BPSK and QPSK
    K0 = 1; % SamplesPerSymbol
    phi = zeros(1,len);
    lambda = zeros(1,len);
    theta = Bn/(xi+1/(4*xi));
    d = 1+2*xi*theta+theta^2;
    g1 = 4*(theta^2/d)/(Kp*K0);
    gp = 4*xi*(theta/d)/(Kp*K0);
    for n  = 2:len
        % phase shift
        y(n) = x(n)*exp(-1i*lambda(n-1));
        % PLL
        err(n) = sign(real(y(n)))*imag(y(n))-sign(imag(y(n)))*real(y(n));
%         err(n) = sign(real(x(n)))*imag(x(n)); % uncomment for BPSK
        phi(n) = g1*err(n) + phi(n-1);
        % DDS
        lambda(n) = (gp*err(n-1)+phi(n-1))+lambda(n-1);
    end
        
end