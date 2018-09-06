%% lfp_dfa
% analyses lfps recorded from several layers across mouse V1
% during periods of locomotion and periods of stationarity.
% Output generated is in 'out' structure. Contains DFA, power and
% amplitudes computed across various frequencies of interest.


% --------------------------
% VERSION 1
% --------------------------
v = 2;
foi_range       = unique(round(2.^[1:.1:9]));
freqs = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
freqs = freqs(1:end-1,:);
dfa_win = [2 100];
% --------------------------

load ~/lfp/dat/expInfo_ori.mat
outdir = '~/lfp/proc/';
run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m

fs = 1250;
siginfo = nbt_Info;
siginfo.converted_sample_frequency = fs;

%% REST

addpath('/home/tpfeffer/Documents/MATLAB/Colormaps/Colormaps (5)/Colormaps/')
load ~/lfp/expInfo.mat
d = dir(sprintf('~/lfp/proc/*dfa*v%d*mat',1));
clear allpow allout

totaldepth = [min([s.series_depth]) max([s.series_depth])];
totaldepth = totaldepth(1)-(25*32):totaldepth(2);

cnt1 = 0; cnt2 = 0;
for idir = 1 : length(d)
  load(['~/lfp/proc/' d(idir).name])
  
  cond(idir)=isempty(strfind(out.mouse_notes,'+/-'));
% 
  if cond(idir) == 1
    cnt1 = cnt1 + 1;
    allout(:,:,1,cnt1) = out.dfa;
    allpow(:,:,1,cnt1) = out.pow;
    allspeed(cnt1,1) = out.mean_speed;
%     alldepth(cnt1,1) = depth;
  else
    cnt2 = cnt2 + 1;
    allout(:,:,2,cnt2) = out.dfa;
    allpow(:,:,2,cnt2) = out.pow;
    allspeed(cnt2,2) = out.mean_speed;
%     alldepth(cnt2,2) = depth;
  end
  
end

close;

dfa_rest = allout;
pow_rest = allpow;

clear allout allpow 

load ~/lfp/expInfo.mat
d = dir(sprintf('~/lfp/proc/*stim*dfa*v%d*mat',2));

cnt1 = 0; cnt2 = 0;
for idir = 1 : length(d)
  load(['~/lfp/proc/' d(idir).name])
  
  cond(idir)=isempty(strfind(out.mouse_notes,'+/-'));
% 
  if cond(idir) == 1
    cnt1 = cnt1 + 1;
    allout(:,:,1,cnt1) = out.dfa;
    allpow(:,:,1,cnt1) = out.pow;
    allspeed(cnt1,1) = out.mean_speed;
  else
    cnt2 = cnt2 + 1;
    allout(:,:,2,cnt2) = out.dfa;
    allpow(:,:,2,cnt2) = out.pow;
    allspeed(cnt2,2) = out.mean_speed;
  end
  
end

close;

out.dfa_stim = allout(:,:,:,1:46);
out.pow_stim = allpow(:,:,:,1:46);

out.dfa_rest = dfa_rest;
out.pow_rest = pow_rest;
%%

% --------------------------
% PLOT CONTRASTS - DFA
% --------------------------

figure; set(gcf,'color','w')

subplot(4,4,[5 6])
par = mean(mean(out.dfa_stim(:,:,1,:),1),4)-mean(mean(out.dfa_rest(:,:,1,:),1),4);
plot(par,'linewidth',3); axis tight; box off
ylabel('DFA exponent'); %xlabel('Frequency')
axis([1 64 -0.1 0.1])
line([1 64],[0 0],'linestyle','--','color','k')
tp_editplots
set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',[])

subplot(2,2,3)

subplot(4,4,[9 10 13 14])
% [~,p] = ttest(out.dfa_stim(:,:,1,:),out.dfa_rest(:,:,1,:),'dim',4);
% [~,p_mask]=fdr(p,0.05);
par = mean(out.dfa_stim(:,:,1,:),4)-mean(out.dfa_rest(:,:,1,:),4);
% par(p_mask==0)=nan;
imagescnan(par,[-0.1 0.1]); colormap(parula)
xlabel('Frequency [Hz]'); ylabel('Channel #')

set(gca,'YTick',1:5:32,'YTickLabels',1:5:32)
set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',mean(out.freqs(1:9:end-1,:),2))
tp_editplots

subplot(4,4,[11 15])
par = flipud(mean(mean(out.dfa_stim(:,:,1,:),4),2))-flipud(mean(mean(out.dfa_rest(:,:,1,:),4),2));
line([0 0],[1 32],'linestyle','--','color','k'); hold on
plot(par,1:32,'linewidth',3); axis tight; box off
xlabel('DFA exponent'); %ylabel('Channel #')
axis([-0.05 0.05 0 32])
set(gca,'YTick',1:5:32,'YTickLabels',[])

tp_editplots

print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_dfa_stimvsrest_v%d.pdf',v))

%%
% --------------------------
% PLOT CONTRASTS - POWER
% --------------------------

figure; set(gcf,'color','w')

subplot(4,4,[5 6])
par = (mean(mean(out.pow_stim(:,:,1,:),1),4)-mean(mean(out.pow_rest(:,:,1,:),1),4))./mean(mean(out.pow_stim(:,:,1,:),1),4);
plot(par,'linewidth',3); axis tight; box off
ylabel('DFA exponent'); %xlabel('Frequency')
axis([1 64 -0.5 1])
line([1 64],[0 0],'linestyle','--','color','k')
tp_editplots
set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',[])

subplot(2,2,3)

subplot(4,4,[9 10 13 14])
% [~,p] = ttest(out.dfa_stim(:,:,1,:),out.dfa_rest(:,:,1,:),'dim',4);
% [~,p_mask]=fdr(p,0.05);
par = (mean(out.pow_rest(:,:,1,:),4)-mean(out.pow_stim(:,:,1,:),4))./mean(out.pow_rest(:,:,1,:),4);
% par(p_mask==0)=nan;
imagescnan(par,[-2.5 5.5]); colormap(parula)
xlabel('Frequency [Hz]'); ylabel('Channel #')

set(gca,'YTick',1:5:32,'YTickLabels',1:5:32)
set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',mean(out.freqs(1:9:end-1,:),2))
tp_editplots

subplot(4,4,[11 15])
par = (flipud(mean(mean(out.pow_stim(:,:,1,:),4),2))-flipud(mean(mean(out.pow_rest(:,:,1,:),4),2)))./flipud(mean(mean(out.pow_stim(:,:,1,:),4),2));
line([0 0],[1 32],'linestyle','--','color','k'); hold on
plot(par,1:32,'linewidth',3); axis tight; box off
xlabel('DFA exponent'); %ylabel('Channel #')
axis([-0.5 1 0 32])
set(gca,'YTick',1:5:32,'YTickLabels',[])

tp_editplots

print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_dfa_stimvsrest_v%d.pdf',v))