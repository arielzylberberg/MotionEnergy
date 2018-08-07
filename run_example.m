
%% example 1: filter plots

load high_coh_right
% rearrange dimensions
stim = permute(S,[3,2,1]);%time,irrel,rel
sf = 75;
motion_energy = motionenergy.motionfilter(stim,sf,'sigmac',0.5,'doplot',1);

%% example 2: a few trials
files = {'high_coh_right','low_coh_right','low_coh_left','high_coh_left'};
sf = 75;
for i=1:length(files)
    load(files{i});
    stim = permute(S,[3,2,1]);%time,irrel,rel
    motion_energy = motionenergy.motionfilter(stim,sf,'sigmac',0.5,'doplot',0);
    ME(:,i) = motion_energy;
end

figure
plot(ME)
legend('R high','R low','L low','L high')
xlabel('time step')
ylabel('motion energy (a.u.)')


