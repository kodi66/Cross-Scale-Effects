% Eissa et al., The Cross-Scale Effects of Neural Interactions during Human 
% Neocortical Seizure Activity, PNAS, 2017.

% load SP
% load rLFP

LFP = zeros(50,50,10000);

order = 3;
fcutlow  = 2;       % [Hz]
fcuthigh = 50;      % [Hz]
fs = 1000;          % sample frequency

[b,a] = butter(order,[fcutlow,fcuthigh]/(fs/2), 'bandpass');

for i = 1:50
    
    for j = 1:50
   
        ts = rLFP(i,j,:);
        LFP(i,j,:) = filtfilt(b,a,ts);
        
    end
    
end

LFPc = zeros(1,10000);      % LFP at the center of the spatial domain
LFPc(:) = LFP(25,25,:);

figure

for t = 1:10000
    
   subplot(3,4,[1 2 5 6])
   imagesc([-2 2],[-2 2],SP(:,:,t))
   shading interp
   axis square
   title('Spiking Activity')
   caxis([0 1])
   
   subplot(3,4,[3 4 7 8])
   imagesc([-2 2],[-2 2],LFP(:,:,t))
   shading interp
   axis square
   title('Synthetic LFP')
   caxis([-150 150])
   
   subplot(3,4,[9 10 11 12])
   box on
   cla
   hold on
   plot(0.001:0.001:10, LFPc)
   plot([t*0.001 t*0.001], [-200 200],'r')
   title('LFP at the center of the spatial domain')
   xlabel('Time [s]')
   drawnow
   
end
