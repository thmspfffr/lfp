%% lfp_dfa
% analyses lfps recorded from several layers across mouse V1
% during periods of locomotion and periods of stationarity.
% Output generated is in 'out' structure. Contains DFA, power and
% amplitudes computed across various frequencies of interest.

clear

% --------------------------
% VERSION 1
% --------------------------
v = 2;
foi_range       = unique(round(2.^[1:.1:9]));
freqs = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
freqs = freqs(1:end-1,:);
dfa_win = [2 100];
% --------------------------

load ~/lfp/expInfo.mat
outdir = '~/lfp/proc/';
run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m

fs = 1250;
siginfo = nbt_Info;
siginfo.converted_sample_frequency = fs;
all_depth = [min([d(:).rel_depth]):25:max([d(:).rel_depth])];

for idir = 1 : length(l)
  
  if ~exist(sprintf([outdir 'lfp_dfa_d%d_v%d_processing.txt'],idir,v))
    system(['touch ' outdir sprintf('lfp_dfa_d%d_v%d_processing.txt',idir,v)]);
  else
    continue
  end
  
  clear out
  out.freqs = freqs;
  
  out.exp_num       = l(idir).exp_num;
  out.mouse_counter = l(idir).mouse_counter;
  out.mouse_id      = m([m.mouse_counter]==out.mouse_counter).mouse_id;
  out.series_num    = l(idir).series_num;
  out.mouse_notes   = m([m.mouse_counter]==out.mouse_counter).mouse_notes;
  out.mean_speed    = l(idir).mean_speed;
  out.mean_speed    = l(idir).percent_moving;
  
  dd = dir(['~/lfp/dat/' sprintf('%s_s%02d_*_%02d.lfp',out.mouse_id, out.series_num,out.exp_num)]);
  
  fid = fopen(['~/lfp/dat/' dd.name]);
  lfpData = fread(fid,'int16');
  fclose(fid);
  
  channel_count = 32;
  % reshape LFP data into channels x samples
  lfpData = reshape(lfpData, channel_count, []);
  
  % relative depths from laura
  rel_depth = [min([d(idir).rel_depth]):25:max([d(idir).rel_depth])];
  indiv_depth_idx = ~isnan(d(idir).rel_depth);
  abs_depth_idx = all_depth>=rel_depth(1) & all_depth<=rel_depth(end);
  
  % preallocate outputs
  out.dfa  = nan(size(all_depth(:),1),length(freqs));
  out.pow  = nan(size(all_depth(:),1),length(freqs));
  out.ampl = nan(size(all_depth(:),1),length(freqs));
  
  for ifoi = 1 : length(freqs)
    
    fprintf('Freq%d...\n',ifoi)
    
    ampenv = single(nbt_filter_fir(lfpData',freqs(ifoi,1),freqs(ifoi,2),siginfo.converted_sample_frequency,2/freqs(ifoi,1)));
    ampenv = abs(hilbert(ampenv));
    
    tmp = tp_dfa(ampenv,dfa_win,fs,0.5,15);
    dfa(:,ifoi) = tmp.exp;
    
    [p,f]=pwelch(lfpData',hanning(4000),2000,[freqs(ifoi,1):0.1:freqs(ifoi,2)],fs);
    
    ampl(:,ifoi)  = mean(ampenv);
    pow(:,ifoi)   = mean(p,1);
    %
  end
  
  out.dfa(abs_depth_idx,:)  = dfa(indiv_depth_idx,:);
  out.pow(abs_depth_idx,:)  = pow(indiv_depth_idx,:);
  out.ampl(abs_depth_idx,:) = ampl(indiv_depth_idx,:);
  
  save(sprintf([outdir 'lfp_dfa_d%d_v%d.mat'],idir,v),'out');
  
end

error('!')

%%
addpath('/home/tpfeffer/Documents/MATLAB/Colormaps/Colormaps (5)/Colormaps/')
load ~/lfp/expInfo.mat
d = dir(sprintf('~/lfp/proc/*dfa*v%d*mat',1));
clear allpow allout

totaldepth = [min([s.series_depth]) max([s.series_depth])];
totaldepth = totaldepth(1)-(25*32):totaldepth(2);
% 
% cnt1 = 0; cnt2 = 0;
% dfa = nan(length(totaldepth),64,124);
% pow = nan(length(totaldepth),64,124);
% 
% for idir = 1 : length(d)
%   idir
%   load(['~/lfp/proc/' d(idir).name])
%   
%   % align based on depth
%   depth = s([s.mouse_counter]==out.mouse_counter & [s.series_num]==out.series_num).series_depth;
%   dd = depth-(25*31):25:depth; 
%   
%   for i = 1 : length(totaldepth)
% %     i
%      curr_depth = [totaldepth(i)-25 totaldepth(i)+25];
%      
%      if curr_depth(1)<min(dd)
%        dfa(i,:,idir) = nan([1 size(out.dfa,2)]);
%        pow(i,:,idir) = nan([1 size(out.pow,2)]);
%      elseif curr_depth(1)>max(dd)
%        dfa(i,:,idir) = nan([1 size(out.dfa,2)]);
%        pow(i,:,idir) = nan([1 size(out.pow,2)]);
%      else
%        idx = dd>=curr_depth(1)&dd<=curr_depth(2);
%        dfa(i,:,idir) = mean(out.dfa(idx,:));
%        pow(i,:,idir) = mean(out.pow(idx,:));
%      end 
%   endwq
% end

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

% cnt1 = 0; cnt2 = 0;
% for idir = 1 : length(d)
%   
%   load(['~/lfp/proc/' d(idir).name])
%   
%   cond(idir)=out.mean_speed<median([l.percent_moving]);
%   
%   if cond(idir) == 1
%     cnt1 = cnt1 + 1;
%      
%     allout(:,:,1,cnt1) = out.dfa;
%     allpow(:,:,1,cnt1) = out.pow;
%   else
%     cnt2 = cnt2 + 1;
%     allout(:,:,2,cnt2) = out.dfa;
%     allpow(:,:,2,cnt2) = out.pow;
%   end
%   
% end


close;

out.dfa = allout;
out.pow = allpow;
%%
% controls

for icond = 1 : 2
  
  figure; set(gcf,'color','w')
  
  subplot(4,4,[5 6])
  par = mean(mean(out.dfa(:,:,icond,:),1),4);
  plot(par,'linewidth',3); axis tight; box off
  ylabel('DFA exponent'); %xlabel('Frequency')
  axis([1 64 min(par)-0.05 max(par)+0.05])
  tp_editplots
  set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',[])
  
  subplot(2,2,3)
  
  subplot(4,4,[9 10 13 14])
  par = mean(out.dfa(:,:,icond,:),4);
  imagesc(par,[0.5 1]); colormap(inferno)
  xlabel('Frequency [Hz]'); ylabel('Channel #')
  
  set(gca,'YTick',1:5:32,'YTickLabels',1:5:32)
  set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',mean(out.freqs(1:9:end-1,:),2))
  tp_editplots
  
  subplot(4,4,[11 15])
  par = flipud(mean(mean(out.dfa(:,:,icond,:),4),2));
  plot(par,1:32,'linewidth',3); axis tight; box off
  xlabel('DFA exponent'); %ylabel('Channel #')
  % axis([min(mean(out.dfa,2))-0.05 max(mean(out.dfa,2))+0.05 0 32])
  set(gca,'YTick',1:5:32,'YTickLabels',[])
  
  tp_editplots
  
  print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_dfa_c%d_v%d.pdf',icond,v))
  
end

% --------------------------
% PLOT CONTRASTS
% --------------------------

figure; set(gcf,'color','w')

subplot(4,4,[5 6])
par = mean(mean(out.dfa(:,:,2,:),1),4)-mean(mean(out.dfa(:,:,1,:),1),4);
plot(par,'linewidth',3); axis tight; box off
ylabel('DFA exponent'); %xlabel('Frequency')
axis([1 64 -0.1 0.1])
line([1 64],[0 0],'linestyle','--','color','k')
tp_editplots
set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',[])

subplot(2,2,3)

subplot(4,4,[9 10 13 14])
[~,p] = ttest(out.dfa(:,:,2,:),out.dfa(:,:,1,:),'dim',4);
[~,p_mask]=fdr(p,0.05);
par = mean(out.dfa(:,:,2,:),4)-mean(out.dfa(:,:,1,:),4);
par(p_mask==0)=nan;
imagescnan(par,[-0.1 0.1]); colormap(parula)
xlabel('Frequency [Hz]'); ylabel('Channel #')

set(gca,'YTick',1:5:32,'YTickLabels',1:5:32)
set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',mean(out.freqs(1:9:end-1,:),2))
tp_editplots

subplot(4,4,[11 15])
par = flipud(mean(mean(out.dfa(:,:,2,:),4),2))-flipud(mean(mean(out.dfa(:,:,1,:),4),2));
line([0 0],[1 32],'linestyle','--','color','k'); hold on
plot(par,1:32,'linewidth',3); axis tight; box off
xlabel('DFA exponent'); %ylabel('Channel #')
axis([-0.05 0.05 0 32])
set(gca,'YTick',1:5:32,'YTickLabels',[])

tp_editplots

print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_dfa_contrast_v%d.pdf',v))

%% POWER

for icond = 1 : 2
  
  figure; set(gcf,'color','w')
  
  subplot(4,4,[5 6])
  par = mean(mean(out.pow(:,:,icond,:),1),4);
  plot(log10(par),'linewidth',3); axis tight; box off
  ylabel('Power (log10)'); %xlabel('Frequency')
%   axis([1 64 min(log10(par)()-0.05 max(log10(par))+0.05])
  tp_editplots
  set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',[])
  
  subplot(2,2,3)
  
  subplot(4,4,[9 10 13 14])
  par = mean(out.pow(:,:,icond,:),4);
  imagesc(log10(par)); colormap(inferno)
  xlabel('Frequency [Hz]'); ylabel('Channel #')
  
  set(gca,'YTick',1:5:32,'YTickLabels',1:5:32)
  set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',mean(out.freqs(1:9:end-1,:),2))
  tp_editplots
  
  subplot(4,4,[11 15])
  par = flipud(mean(mean(out.pow(:,:,icond,:),4),2));
  plot(par,1:32,'linewidth',3); axis tight; box off
  xlabel('Power (log10)'); %ylabel('Channel #')
  % axis([min(mean(out.dfa,2))-0.05 max(mean(out.dfa,2))+0.05 0 32])
  set(gca,'YTick',1:5:32,'YTickLabels',[])
  
  tp_editplots
  
  print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_pow_c%d_v%d.pdf',icond,v))
  
end

% --------------------------
% PLOT CONTRASTS
% --------------------------

figure; set(gcf,'color','w')

subplot(4,4,[5 6])
par = mean(mean(out.pow(:,:,2,:),1),4)./mean(mean(out.pow(:,:,1,:),1),4);
plot(par,'linewidth',3); axis tight; box off
ylabel('Diff. (Power)'); %xlabel('Frequency')
line([1 64],[1 1],'linestyle','--','color','k')
axis([1 64 0.5 1.5])
tp_editplots
set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',[])

subplot(2,2,3)

subplot(4,4,[9 10 13 14])
[~,p] = ttest(out.pow(:,:,2,:),out.pow(:,:,1,:),'dim',4);
[~,p_mask]=fdr(p,0.05);
par = mean(out.pow(:,:,2,:),4)./mean(out.pow(:,:,1,:),4);
par(p_mask==0)=nan;
imagescnan(par,[0.5 1.5]); colormap(parula)
xlabel('Frequency [Hz]'); ylabel('Channel #')

set(gca,'YTick',1:5:32,'YTickLabels',1:5:32)
set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',mean(out.freqs(1:9:end-1,:),2))
tp_editplots

subplot(4,4,[11 15])
par = flipud(mean(mean(out.pow(:,:,2,:),4),2))./flipud(mean(mean(out.pow(:,:,1,:),4),2));
plot(par,1:32,'linewidth',3); axis tight; box off
xlabel('Diff. (Power)'); %ylabel('Channel #')
axis([0.5 1.5 0 32])
set(gca,'YTick',1:5:32,'YTickLabels',[])
line([1 1],[1 64],'linestyle','--','color','k')

tp_editplots

print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_pow_contrast_v%d.pdf',v))

%% POWER DFA
for icond = 1 : 2
  for ii = 1 : 62

    tmp = out.dfa(:,:,icond,ii);
    dfa(:,ii,icond) = tmp(p_mask);
    tmp = out.pow(:,:,icond,ii);
    pow(:,ii,icond) = tmp(p_mask);
  
  end
end

