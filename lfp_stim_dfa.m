%% lfp_stim_dfa
% analyses lfps recorded from several layers across mouse V1
% during periods of locomotion and periods of stationarity.
% Output generated is in 'out' structure. Contains DFA, power and
% amplitudes computed across various frequencies of interest.

clear 
% --------------------------
% VERSION 1
% --------------------------
v = 1;
foi_range       = unique(round(2.^[1:.1:8]));
freqs = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
freqs = freqs(1:end-1,:);
dfa_win = [2 100];
trl_based= 0;
% --------------------------

load ~/lfp/dat/expInfo_ori.mat
outdir = '~/lfp/proc/';
run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m

fs = 1250;
siginfo = nbt_Info;
siginfo.converted_sample_frequency = fs;
%%
% figure;
warning('off')
for iexp = 1 : length(e)
  
  if ~exist(sprintf([outdir 'lfp_stim_dfa_exp%d_v%d_processing.txt'],iexp,v))
    system(['touch ' outdir sprintf('lfp_stim_dfa_exp%d_v%d_processing.txt',iexp,v)]);
  else
    continue
  end
  
  out.exp_num       = e(iexp).exp_num;
  out.series_num    = e(iexp).series_num;
  out.mouse_counter = e(iexp).mouse_counter;
  out.mouse_id      = m([m.mouse_counter]==out.mouse_counter).mouse_id;
  out.mouse_notes   = m([m.mouse_counter]==out.mouse_counter).mouse_notes;
  out.cond          = 2-isempty(strfind(out.mouse_notes,'+/-'));
  out.cond_descr    = '1=control; 2=pv+knockout';
  
  trls = find([t.exp_num] == out.exp_num & [t.series_num] ==out.series_num & [t.mouse_counter]==out.mouse_counter);
  out.trials = trls;
  
  % load and reshape data
  d = dir(['~/lfp/dat/' sprintf('%s_s%02d_*_%02d.lfp',out.mouse_id, out.series_num,out.exp_num)]);
  fid = fopen(['~/lfp/dat/' d.name]);
  lfpData = fread(fid,'int16');
  fclose(fid);
  lfpData = reshape(lfpData, 32, []);
  
  dat = [];
  
  if trl_based
    dat = zeros(32,2001,length(trls));
    for itrl = 1 : length(trls)
      
      out.time          = [t(trls(itrl)).trial_onset/(30000/fs) t(trls(itrl)).trial_offset/(30000/fs)];
      out.contrast      = c(t(trls(itrl)).trial_num).grat_contrast;      
      epoch_ts = double([int64(out.time(1)-0.3*fs) : int64(out.time(1) + 1.3*fs)]);
      dat(:,:,itrl) = lfpData(:,epoch_ts);

    end
  else
    dat = lfpData; clear lfpData
  end
  
  for ifoi = 1 : length(freqs)
    
    fprintf('Freq%d...\n',ifoi)
    
    ampenv = single(nbt_filter_fir(dat',freqs(ifoi,1),freqs(ifoi,2),siginfo.converted_sample_frequency,2/freqs(ifoi,1)));
    ampenv = abs(hilbert(ampenv));
    
    tmp = tp_dfa(ampenv,dfa_win,fs,0.5,15);
    out.dfa(:,ifoi) = tmp.exp;
    
    [p,f]=pwelch(dat',hanning(4000),2000,[freqs(ifoi,1):0.1:freqs(ifoi,2)],fs);
    
    out.ampl(:,ifoi)  = mean(ampenv);
    out.pow(:,ifoi)   = mean(p,1);
    
    % Compute autocorrelations
%     for i = 1 : 32
%       [acf,lags]=autocorr(ampenv(:,i),fs);
%       out.lambda(i,ifoi)=1./tp_fitexpdecay(acf,1:fs+1,0.02);
%     end
    %
  end
 
  save(sprintf([outdir 'lfp_stim_dfa_exp%d_v%d.mat'],iexp,v),'out');

end

error('!')


%%
clear out allpow alldfa cond d 

addpath('/home/tpfeffer/Documents/MATLAB/Colormaps/Colormaps (5)/Colormaps/')
load ~/lfp/dat/expInfo_ori.mat
dd = dir(sprintf('~/lfp/proc/lfp_stim_dfa_exp*v%d.mat',1));
 
all_depth = [min([d(:).rel_depth]):25:max([d(:).rel_depth])];

% relative depths from laura
alldfa   = nan(40,size(freqs,1),2,50);
allpow   = nan(40,size(freqs,1),2,50);
allampl  = nan(40,size(freqs,1),2,50);

cnt1 = 0; cnt2 = 0; subj = 1;

for idir = 1 : length(dd)

  load(['~/lfp/proc/' dd(idir).name])
  clear pow
  
%   for ifreq = 1 : length(freqs)
%     pow(:,ifreq) = mean(out.pow(find(out.freq>=freqs(ifreq,1) & out.freq<=freqs(ifreq,2)),:),1);
%   end
  
  % ALIGN RECORDINGS
  depth_idx       = [d.mouse_counter]==out.mouse_counter&[d.series_num]==out.series_num;
  rel_depth       = min([d(depth_idx).rel_depth]):25:max([d(depth_idx).rel_depth]);
  indiv_depth_idx = ~isnan(d(depth_idx).rel_depth);
  abs_depth_idx   = all_depth>=rel_depth(1) & all_depth<=rel_depth(end);

  cond(idir)=isempty(strfind(out.mouse_notes,'+/-'));
    
  if cond(idir) == 1
    cnt1 = cnt1 + 1;
    alldfa(abs_depth_idx,:,1,cnt1)  = out.dfa(indiv_depth_idx,:);
    allpow(abs_depth_idx,:,1,cnt1)  = out.pow(indiv_depth_idx,:);
%     allampl(abs_depth_idx,:,1,cnt1) = out.ampl(indiv_depth_idx,:);
  else
    cnt2 = cnt2 + 1;
    alldfa(abs_depth_idx,:,2,cnt2)  = out.dfa(indiv_depth_idx,:);
    allpow(abs_depth_idx,:,2,cnt2)  = out.pow(indiv_depth_idx,:);
%     allampl(abs_depth_idx,:,2,cnt2) = out.ampl(indiv_depth_idx,:);
  end 
end

out.freqs = freqs;

out.dfa  = alldfa(~isnan(nanmean(alldfa(:,1,1,:),4)) & ~isnan(nanmean(alldfa(:,1,2,:),4)),:,:,:);
out.pow  = allpow(~isnan(nanmean(alldfa(:,1,1,:),4)) & ~isnan(nanmean(alldfa(:,1,2,:),4)),:,:,:);
% out.ampl = allampl(~isnan(nanmean(alldfa(:,1,1,:),4)) & ~isnan(nanmean(alldfa(:,1,2,:),4)),:,:,:);

out.dfa  = out.dfa(:,:,:,1:min([cnt1 cnt2]));
out.pow  = out.pow(:,:,:,1:min([cnt1 cnt2]));
% out.ampl = allampl;

all_depth = all_depth(~isnan(nanmean(alldfa(:,1,1,:),4)) & ~isnan(nanmean(alldfa(:,1,2,:),4)));

%%
load ~/lfp/expInfo.mat

for icond = 1 : 2
  
  figure; set(gcf,'color','w')
  
  subplot(4,4,[5 6])
  par = mean(nanmean(out.dfa(:,:,icond,:),1),4);
  plot(par,'linewidth',3); axis tight; box off
  line([1 54],[mean(mean(nanmean(out.dfa(:,:,icond,:),1),4)) mean(mean(nanmean(out.dfa(:,:,icond,:),1),4))],'linestyle',':')
  ylabel('DFA exponent'); %xlabel('Frequency')
  axis([1 54 0.5 1])
  tp_editplots
  set(gca,'XTick',1:4:size(out.freqs,1),'XTickLabels',[])
  set(gca,'YTick',[0.5 :0.1: 1],'YTickLabels',num2cell([0.5 :0.1: 1.0]))
      
  subplot(4,4,7); hold on; colormap(gca, 'plasma'); c = colorbar; axis off
  c.Position = [0.58 0.55 0.02 0.15]; 
  c.Ticks = [0 1];
  clim = ['0.5';'1.0'];
  c.TickLabels = clim;
  
  subplot(4,4,[9 10 13 14]); hold on
  par = nanmean(out.dfa(:,:,icond,:),4);
  imagesc(par,[str2num(clim(1,:)) str2num(clim(2,:))]); colormap(plasma)
  xlabel('Frequency [Hz]'); ylabel('Channel #')
  axis([1 54 1 37])
  set(gca,'YTick',1:4:length(all_depth),'YTickLabels',num2cell(all_depth(1:4:end)))
  set(gca,'XTick',1:4:size(out.freqs,1),'XTickLabels',mean(out.freqs(1:4:end,:),2))
  % layers
  line([1 54],[find(all_depth==-60) find(all_depth==-60)],'linestyle','--','color','w')
  line([1 54],[find(all_depth==65) find(all_depth==65)],'linestyle','--','color','w')
  tp_editplots
  set(gca,'YDir','reverse')
  
  subplot(4,4,[11 15])
  par = flipud(mean(nanmean(out.dfa(:,:,icond,:),4),2));
  plot(par,1:size(out.dfa,1),'linewidth',3); axis tight; box off
  xlabel('DFA exponent'); 
  set(gca,'YTick',1:4:length(all_depth),'YTickLabels',[])
  tp_editplots
  axis([0.5 1 1 size(out.dfa,1)])

  print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_dfa_stim_c%d_v%d.pdf',icond,v))
  
end

% --------------------------
% PLOT CONTRASTS
% --------------------------

figure; set(gcf,'color','w')

subplot(4,4,[5 6])
par = mean(nanmean(out.dfa(:,:,2,:),1),4)-mean(nanmean(out.dfa(:,:,1,:),1),4);
plot(par,'linewidth',3); axis tight; box off
ylabel('DFA exponent'); %xlabel('Frequency')
axis([1 54 -0.2 0.2])
line([1 54],[0 0],'linestyle',':','color','k')
tp_editplots
set(gca,'XTick',1:4:size(out.freqs,1),'XTickLabels',[])

subplot(4,4,7); hold on; colormap(gca, 'parula'); c = colorbar; axis off
c.Position = [0.58 0.55 0.02 0.15];
c.Ticks = [0 1];
clim = ['0.1';'0.1'];
c.TickLabels = clim;

subplot(4,4,[9 10 13 14]); hold on
[~,p] = ttest(out.dfa(:,:,2,:),out.dfa(:,:,1,:),'dim',4);
[~,p_mask]=fdr(p,0.05);
par = nanmean(out.dfa(:,:,2,:),4)-nanmean(out.dfa(:,:,1,:),4);
% par(p_mask==0)=nan;
imagescnan(par,[-str2num(clim(1,:)) str2num(clim(2,:))]); %colormap(parula)
% layers
line([1 54],[find(all_depth==-60) find(all_depth==-60)],'linestyle','--','color','w')
line([1 54],[find(all_depth==65) find(all_depth==65)],'linestyle','--','color','w')
xlabel('Frequency [Hz]'); ylabel('Channel #')
set(gca,'YTick',1:4:length(all_depth),'YTickLabels',num2cell(all_depth(1:4:end)))
set(gca,'XTick',1:4:size(out.freqs,1),'XTickLabels',mean(out.freqs(1:4:end,:),2))
axis([1 54 1 37])
tp_editplots
set(gca,'YDir','reverse')

subplot(4,4,[11 15])
par = flipud(mean(nanmean(out.dfa(:,:,2,:),4),2))-flipud(mean(nanmean(out.dfa(:,:,1,:),4),2));
line([0 0],[1 size(out.dfa,1)],'linestyle','--','color','k'); hold on
plot(par,1:size(out.dfa,1),'linewidth',3); axis tight; box off
xlabel('DFA exponent'); %ylabel('Channel #')
set(gca,'YTick',1:4:length(all_depth),'YTickLabels',[])
axis([-0.2 0.2 0 37])

tp_editplots

print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_dfa_stim_contrast_v%d.pdf',v))

%% POWER 

for icond = 1 : 2
  
  figure; set(gcf,'color','w')
  
  subplot(4,4,[5 6]); hold on
  par = mean(nanmean(out.pow(:,:,icond,:),1),4);
  plot(log10(par),'linewidth',3); axis tight; box off
  ylabel('Power (log10)'); %xlabel('Frequency')
  set(gca,'XTick',1:4:size(out.freqs,1),'XTickLabels',[])
  set(gca,'YTick',1:5,'YTickLabels',1:5)
  tp_editplots
  axis([1 size(out.pow,2) 1 5])
  
  subplot(4,4,7); hold on; colormap(gca, 'plasma'); c = colorbar; axis off
  c.Position = [0.58 0.55 0.02 0.15]; 
  c.Ticks = [0 1];
  clim = ['1';'5'];
  c.TickLabels = clim;
  
  subplot(4,4,[9 10 13 14]); hold on
  par = nanmean(out.pow(:,:,icond,:),4);
  imagesc(log10(par),[str2num(clim(1)) str2num(clim(2))]); %colormap(inferno)
  xlabel('Frequency [Hz]'); ylabel('Channel #')
  % layers
  line([1 size(out.pow,2)],[find(all_depth==-60) find(all_depth==-60)],'linestyle','--','color','w')
  line([1 size(out.pow,2)],[find(all_depth==65) find(all_depth==65)],'linestyle','--','color','w')
  set(gca,'YTick',1:4:length(all_depth),'YTickLabels',num2cell(all_depth(1:4:end)))
  set(gca,'XTick',1:4:size(out.freqs,1),'XTickLabels',num2cell(mean(out.freqs(1:4:end,:),2)))
  axis([1 size(out.pow,2) 1 size(out.dfa,1)])
  tp_editplots; set(gca,'YDir','reverse')
  colormap(plasma)
  
  subplot(4,4,[11 15])
  par = flipud(mean(nanmean(out.pow(:,:,icond,:),4),2));
  plot(log10(par),1:size(out.pow,1),'linewidth',3); axis tight; box off
  xlabel('Power (log10)'); %ylabel('Channel #')
  set(gca,'YTick',1:4:length(all_depth),'YTickLabels',[])
  set(gca,'XTick',1:5,'XTickLabels',1:5)
  axis([1 5 1 size(out.pow,1)])
  
  tp_editplots
  
  print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_pow_stim_c%d_v%d.pdf',icond,v))
  
end

% --------------------------
% PLOT CONTRASTS
% --------------------------

figure; set(gcf,'color','w')

subplot(4,4,[5 6]); hold on
par = 100*(mean(nanmean(out.pow(:,:,2,:),1),4)-mean(nanmean(out.pow(:,:,1,:),1),4))./mean(nanmean(out.pow(:,:,1,:),1),4);
plot(par,'linewidth',3); axis tight; box off
ylabel('Change in power'); %xlabel('Frequency')
line([1 size(out.pow,2)],[1 1],'linestyle',':','color','k')
axis([1 size(out.pow,2) -100 100])
tp_editplots
set(gca,'XTick',1:4:size(out.freqs,1),'XTickLabels',[])

subplot(4,4,7); hold on; colormap(gca, 'parula'); c = colorbar; axis off
c.Position = [0.58 0.55 0.02 0.15];
c.Ticks = [0 1];
clim = ['100';'100'];
c.TickLabels = clim;

subplot(4,4,[9 10 13 14]); hold on
% [~,p] = ttest(out.pow(:,:,2,:),out.pow(:,:,1,:),'dim',4);
% [~,p_mask]=fdr(p,0.05);
par = 100*(nanmean(out.pow(:,:,2,:),4)- nanmean(out.pow(:,:,1,:),4))./nanmean(out.pow(:,:,1,:),4);
% par(p_mask==0)=nan;
imagescnan(par,[-str2num(clim(1,:)) str2num(clim(2,:))]); %colormap(redblue)
% layers
line([1 size(out.pow,2)],[find(all_depth==-60) find(all_depth==-60)],'linestyle','--','color','w')
line([1 size(out.pow,2)],[find(all_depth==65) find(all_depth==65)],'linestyle','--','color','w')
xlabel('Frequency [Hz]'); ylabel('Channel #')
set(gca,'YDir','reverse')
axis([1 size(out.pow,2) 1 size(out.dfa,1)])
set(gca,'YTick',1:4:length(all_depth),'YTickLabels',num2cell(all_depth(1:4:end)))
set(gca,'XTick',1:4:size(out.freqs,1),'XTickLabels',mean(out.freqs(1:4:end,:),2))
tp_editplots
% colormap(parula)

subplot(4,4,[11 15]); hold on
par = 100*(flipud(mean(nanmean(out.pow(:,:,2,:),4),2))-flipud(mean(nanmean(out.pow(:,:,1,:),4),2)))./flipud(mean(nanmean(out.pow(:,:,1,:),4),2));
plot(par,1:size(out.pow,1),'linewidth',3); axis tight; box off
xlabel('Change in power');
set(gca,'YTick',1:4:size(out.pow,1),'YTickLabels',[])
set(gca,'XTick',[-100 0 100],'XTickLabels',[-100 0 100])
axis([-100 100 1 size(out.pow,1)])

line([1 1],[1 size(out.pow,1)],'linestyle',':','color','k')

tp_editplots

print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_pow_stim_contrast_v%d.pdf',v))



%% AMPLITUDE

for icond = 1 : 2
  
  figure; set(gcf,'color','w')
  
  subplot(4,4,[5 6])
  par = mean(nanmean(out.ampl(:,:,icond,:),1),4);
  plot(log10(par),'linewidth',3); axis tight; box off
  ylabel('Power (log10)'); %xlabel('Frequency')
%   axis([1 64 min(log10(par)()-0.05 max(log10(par))+0.05])
  tp_editplots
  set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',[])

  subplot(2,2,3)
  
  subplot(4,4,[9 10 13 14])
  par = nanmean(out.ampl(:,:,icond,:),4);
  imagesc(log10(par)); colormap(plasma)
  xlabel('Frequency [Hz]'); ylabel('Channel #')
  
  set(gca,'YTick',1:4:length(all_depth),'YTickLabels',num2cell(all_depth(1:4:end)))
  set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',mean(out.freqs(1:9:end-1,:),2))
  tp_editplots
  
  subplot(4,4,[11 15])
  par = flipud(mean(nanmean(out.ampl(:,:,icond,:),4),2));
  plot(par,1:size(out.ampl,1),'linewidth',3); axis tight; box off
  xlabel('Amplitude (log10)'); %ylabel('Channel #')
  % axis([min(mean(out.dfa,2))-0.05 max(mean(out.dfa,2))+0.05 0 32])
  set(gca,'YTick',1:4:length(all_depth),'YTickLabels',[])
  
  tp_editplots
  
  print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_pow_c%d_v%d.pdf',icond,v))
  
end

% --------------------------
% PLOT CONTRASTS
% --------------------------

figure; set(gcf,'color','w')

subplot(4,4,[5 6])
par = mean(nanmean(out.ampl(:,:,2,:),1),4)./mean(nanmean(out.ampl(:,:,1,:),1),4);
plot(par,'linewidth',3); axis tight; box off
ylabel('Diff. (Power)'); %xlabel('Frequency')
line([1 54],[1 1],'linestyle','--','color','k')
axis([1 54 0.66 1.5])
tp_editplots
set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',[])

subplot(2,2,3)

subplot(4,4,[9 10 13 14])
% [~,p] = ttest(out.pow(:,:,2,:),out.pow(:,:,1,:),'dim',4);
% [~,p_mask]=fdr(p,0.05);
par = nanmean(out.ampl(:,:,2,:),4)./nanmean(out.ampl(:,:,1,:),4);
% par(p_mask==0)=nan;
imagescnan(par,[0.66 1.5]); colormap(redblue)
xlabel('Frequency [Hz]'); ylabel('Channel #')

set(gca,'YTick',1:4:length(all_depth),'YTickLabels',num2cell(all_depth(1:4:end)))
set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',mean(out.freqs(1:9:end-1,:),2))
tp_editplots

subplot(4,4,[11 15])
par = flipud(mean(nanmean(out.ampl(:,:,2,:),4),2))./flipud(mean(nanmean(out.ampl(:,:,1,:),4),2));
plot(par,1:size(out.ampl,1),'linewidth',3); axis tight; box off
xlabel('Diff. (Power)'); %ylabel('Channel #')
axis([0.1 2 0 size(out.ampl,1)])
set(gca,'YTick',1:5:40,'YTickLabels',[])
line([1 1],[1 64],'linestyle','--','color','k')

tp_editplots

print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_ampl_contrast_v%d.pdf',v))

%% POWER DFA
for icond = 1 : 2
  for ii = 1 : 62

    tmp = out.dfa(:,:,icond,ii);
    dfa(:,ii,icond) = tmp(p_mask);
    tmp = out.pow(:,:,icond,ii);
    pow(:,ii,icond) = tmp(p_mask);
  
  end
end
