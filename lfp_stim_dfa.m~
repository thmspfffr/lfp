%% lfp_stim_dfa
% analyses lfps recorded from several layers across mouse V1
% during periods of locomotion and periods of stationarity.
% Output generated is in 'out' structure. Contains DFA, power and
% amplitudes computed across various frequencies of interest.


% --------------------------
% VERSION 2
% --------------------------
v = 2;
foi_range 	= unique(round(2.^[1:.1:7.5]));
freqs       = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
freqs       = freqs(1:end-1,:);
dfa_win     = [3 75];
trl_based   = 0;
% --------------------------
% VERSION 3
% --------------------------
v = 3;
foi_range 	= unique(round(2.^[1:.1:7.5]));
freqs       = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
freqs       = freqs(1:end-1,:);
dfa_win     = [3 75];
trl_based   = 1;
% --------------------------

load ~/lfp/dat/expInfo_ori.mat
outdir = '~/lfp/proc/';
run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m

fs = 1250;
siginfo = nbt_Info;
siginfo.converted_sample_frequency = fs;
%%
figure;
warning('off')
for iexp = 1 : length(e)
  
%   if ~exist(sprintf([outdir 'lfp_stim_dfa_exp%d_v%d_processing.txt'],iexp,v))
%     system(['touch ' outdir sprintf('lfp_stim_dfa_exp%d_v%d_processing.txt',iexp,v)]);
%   else
%     continue
%   end
%   
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
% itrl
      out.time          = [t(trls(itrl)).trial_onset/(30000/fs) t(trls(itrl)).trial_offset/(30000/fs)];
      out.contrast      = c(t(trls(itrl)).trial_num).grat_contrast;
%       cnt_contrast(
      
      epoch_ts = double([int64(out.time(1)-0.3*fs) : int64(out.time(1) + 1.3*fs)]);
      
      dat(:,:,itrl) = lfpData(:,epoch_ts);


    end
  else
    dat = lfpData; clear lfpData
  end
  
  plot(mean(mean(dat,3)),'color',[0.8 0.8 0.8]-rand/5); hold on; drawnow
%   resample_fs = fs/2;
%   
%   dat = resample(dat',1,fs/resample_fs)';
%   
%   [out.pow,out.freq]=pwelch(dat',hanning(resample_fs),resample_fs/5,[1:0.2:150],resample_fs);
%   
%   for ifoi = 1 : length(freqs)
%     ifoi
%     ampenv = nbt_filter_fir(dat',freqs(ifoi,1),freqs(ifoi,2),resample_fs,2/freqs(ifoi,1));
%     ampenv = abs(hilbert(ampenv));
%     
%     tmp = tp_dfa(ampenv,dfa_win,resample_fs,0.5,15);
%     
%     out.dfa(:,ifoi) = tmp.exp;
%     for i = 1 : 32
%       [acf,lags]=autocorr(ampenv(:,i),resample_fs);
%       out.lambda(i,ifoi)=1./tp_fitexpdecay(acf,1:resample_fs+1,0.02);
%     end
%   end
%   
%   save(sprintf([outdir 'lfp_stim_dfa_exp%d_v%d.mat'],iexp,v),'out');

end

error('!')

set(gca,'XTick',[125:250:16*125],'XTickLabels',num2cell([round((-0.2:0.2:1.4)*10)/10]))
line([375 375], [-600 800],'color','k','linestyle',':')
line([1625 1625], [-600 800],'color','k','linestyle',':')

%%
iexp =4
load(sprintf([outdir 'lfp_stim_dfa_exp%d_v%d.mat'],iexp,2),'out');

f = out.freq>15&out.freq<50;
for ichan = 1 : 32
  beta(ichan,:)=tp_fitpowlaw(out.pow(f,ichan),out.freq(f),[out.pow(find(f,1,'first'),ichan) 2]);
end
% beta(2)

figure;
plot(beta(:,2),1:32,'.','markersize',20)
set(gca,'YDir','reverse')
axis([-3  1 0 33])
%%
addpath('/home/tpfeffer/Documents/MATLAB/Colormaps/Colormaps (5)/Colormaps/')
load ~/lfp/expInfo.mat

clear out

d = dir(sprintf('~/lfp/proc/*stim*dfa*v%d*mat',2));
d_res = dir(sprintf('~/lfp/proc/lfp_dfa*v%d*mat',2));

clear allpow allout

% totaldepth = [min([s.series_depth]) max([s.series_depth])];
% totaldepth = totaldepth(1)-(25*32):totaldepth(2);
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

clear allpow1 allpow2 allout1 allout2
cnt1 = 0; cnt2 = 0;
for idir = 1 : length(d)
  idir
  load(['~/lfp/proc/' d(idir).name])
  
  cond(idir)=isempty(strfind(out.mouse_notes,'+/-'));
  %
  if cond(idir) == 1
    cnt1 = cnt1 + 1;
    allout1(:,:,cnt1) = out.dfa;
    allpow1(:,:,cnt1) = out.pow;
    alllambda1(:,:,cnt1) = out.lambda;
%     allspeed(cnt1,1) = out.mean_speed;
    %     alldepth(cnt1,1) = depth;
  else
    cnt2 = cnt2 + 1;
    allout2(:,:,cnt2) = out.dfa;
    allpow2(:,:,cnt2) = out.pow;
    alllambda2(:,:,cnt2) = out.lambda;
%     allspeed(cnt2,2) = out.mean_speed;
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

clear out

out.dfa(:,:,:,1) = allout1(:,:,1:40);
out.dfa(:,:,:,2) = allout2(:,:,1:40);
out.lambda(:,:,:,1) = alllambda1(:,:,1:40);

out.pow(:,:,:,1) = allpow1(:,:,1:40);
out.pow(:,:,:,2) = allpow2(:,:,1:40);
out.lambda(:,:,:,2) = alllambda2(:,:,1:40);

%%
% controls
out.freqs = freqs;

for icond = 1 : 2
  
  figure; set(gcf,'color','w')
  
  subplot(4,4,[5 6])
  par = mean(mean(out.dfa(:,:,:,icond),1),3);
  plot(par,'linewidth',3,'color','k'); axis tight; box off
  ylabel('DFA exponent'); %xlabel('Frequency')
  axis([1 64 min(par)-0.05 max(par)+0.05])
  tp_editplots
  set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',[])
  
  subplot(2,2,3)
  
  subplot(4,4,[9 10 13 14])
  par = mean(out.dfa(:,:,:,icond),3);
  imagesc(par,[0.5 1]); colormap(inferno)
  xlabel('Frequency [Hz]'); ylabel('Channel #')
  
  set(gca,'YTick',1:5:32,'YTickLabels',1:5:32)
  set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',mean(out.freqs(1:9:end-1,:),2))
  tp_editplots
  
  subplot(4,4,[11 15])
  par = flipud(mean(mean(out.dfa(:,:,:,icond),3),2));
  plot(par,1:32,'linewidth',3,'color','k'); axis tight; box off
  xlabel('DFA exponent'); %ylabel('Channel #')
  % axis([min(mean(out.dfa,2))-0.05 max(mean(out.dfa,2))+0.05 0 32])
  set(gca,'YTick',1:5:32,'YTickLabels',[])
  
  tp_editplots
  
  print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_stim_dfa_c%d_v%d.pdf',icond,v))
  
end

% --------------------------
% PLOT CONTRASTS
% --------------------------

figure; set(gcf,'color','w')

subplot(4,4,[5 6])
par = mean(mean(out.dfa(:,:,:,2),1),3)-mean(mean(out.dfa(:,:,:,1),1),3);
plot(par,'linewidth',3,'color',[0.7 0.7 0.7]); axis tight; box off
ylabel('DFA exponent'); %xlabel('Frequency')
axis([1 64 -0.15 0.15])
line([1 64],[0 0],'linestyle','--','color','k')
tp_editplots
set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',[])

subplot(2,2,3)

subplot(4,4,[9 10 13 14])
[~,p] = ttest(out.dfa(:,:,:,2),out.dfa(:,:,:,1),'dim',3);
[~,p_mask]=fdr(p,0.05);
par = mean(out.dfa(:,:,:,2),3)-mean(out.dfa(:,:,:,1),3);
par(p_mask==0)=nan;
imagescnan(par,[-0.15 0.15]); colormap(cmap)
xlabel('Frequency [Hz]'); ylabel('Channel #')

set(gca,'YTick',1:5:32,'YTickLabels',1:5:32)
set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',mean(out.freqs(1:9:end-1,:),2))
tp_editplots

subplot(4,4,[11 15])
par = flipud(mean(mean(out.dfa(:,:,:,2),3),2))-flipud(mean(mean(out.dfa(:,:,:,1),3),2));
line([0 0],[1 32],'linestyle','--','color','k'); hold on
plot(par,1:32,'linewidth',3); axis tight; box off
xlabel('DFA exponent'); %ylabel('Channel #')
axis([-0.05 0.05 0 32])
set(gca,'YTick',1:5:32,'YTickLabels',[])

tp_editplots

print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_stim_dfa_contrast_v%d.pdf',v))

%% POWER

for icond = 1 : 2
  
  figure; set(gcf,'color','w')
  
  subplot(4,4,[5 6])
  par = mean(mean(out.pow(:,:,:,icond),1),3);
  plot(log10(par),'linewidth',3); axis tight; box off
  ylabel('Power (log10)'); %xlabel('Frequency')
  %   axis([1 64 min(log10(par)()-0.05 max(log10(par))+0.05])
  tp_editplots
  set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',[])
  
  subplot(2,2,3)
  
  subplot(4,4,[9 10 13 14])
  par = mean(out.pow(:,:,:,icond),3);
  imagesc(log10(par)); colormap(inferno)
  xlabel('Frequency [Hz]'); ylabel('Channel #')
  
  set(gca,'YTick',1:5:32,'YTickLabels',1:5:32)
  set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',mean(out.freqs(1:9:end-1,:),2))
  tp_editplots
  
  subplot(4,4,[11 15])
  par = flipud(mean(mean(out.pow(:,:,:,icond),3),2));
  plot(par,1:32,'linewidth',3); axis tight; box off
  xlabel('Power (log10)'); %ylabel('Channel #')
  % axis([min(mean(out.dfa,2))-0.05 max(mean(out.dfa,2))+0.05 0 32])
  set(gca,'YTick',1:5:32,'YTickLabels',[])
  
  tp_editplots
  
  print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_stim_pow_c%d_v%d.pdf',icond,v))
  
end

% --------------------------
% PLOT CONTRASTS
% --------------------------

figure; set(gcf,'color','w')

subplot(4,4,[5 6])
par = mean(mean(out.pow(:,:,:,2),1),3)./mean(mean(out.pow(:,:,:,1),1),3);
plot(par,'linewidth',3); axis tight; box off
ylabel('Diff. (Power)'); %xlabel('Frequency')
line([1 64],[1 1],'linestyle','--','color','k')
axis([1 64 0.2 1.8])
tp_editplots
set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',[])

subplot(2,2,3)

subplot(4,4,[9 10 13 14])
[~,p] = ttest(out.pow(:,:,:,2),out.pow(:,:,:,1),'dim',3);
[~,p_mask]=fdr(p,0.05);
par = mean(out.pow(:,:,:,2),3)./mean(out.pow(:,:,:,1),3);
par(p_mask==0)=nan;
imagescnan(par,[0.2 1.8]); colormap(cmap)
xlabel('Frequency [Hz]'); ylabel('Channel #')

set(gca,'YTick',1:5:32,'YTickLabels',1:5:32)
set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',mean(out.freqs(1:9:end-1,:),2))
tp_editplots

subplot(4,4,[11 15])
par = flipud(mean(mean(out.pow(:,:,:,2),3),2))./flipud(mean(mean(out.pow(:,:,:,1),3),2));
plot(par,1:32,'linewidth',3); axis tight; box off
xlabel('Diff. (Power)'); %ylabel('Channel #')
axis([0.2 1.8 0 32])
set(gca,'YTick',1:5:32,'YTickLabels',[])
line([1 1],[1 64],'linestyle','--','color','k')

tp_editplots

print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_stim_pow_contrast_v%d.pdf',v))

%% POWER DFA
for icond = 1 : 2
  for ii = 1 : 62
    
    tmp = out.dfa(:,:,icond,ii);
    dfa(:,ii,icond) = tmp(p_mask);
    tmp = out.pow(:,:,icond,ii);
    pow(:,ii,icond) = tmp(p_mask);
    
  end
end

%% CONTRAST

addpath('/home/tpfeffer/Documents/MATLAB/Colormaps/Colormaps (5)/Colormaps/')
load ~/lfp/expInfo.mat

outa=out;

d = dir(sprintf('~/lfp/proc/*stim*dfa*v%d*mat',2));
d_res = dir(sprintf('~/lfp/proc/lfp_dfa*v%d*mat',2));

clear allout1 allout2 allpow2 allpow1

% totaldepth = [min([s.series_depth]) max([s.series_depth])];
% totaldepth = totaldepth(1)-(25*32):totaldepth(2);
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

clear allpow1 allpow2 allout1 allout2
cnt1 = 0; cnt2 = 0;
for idir = 1 : length(d_res)
  idir
  load(['~/lfp/proc/' d_res(idir).name])
  
  cond(idir)=isempty(strfind(out.mouse_notes,'+/-'));
  %
  if cond(idir) == 1
    cnt1 = cnt1 + 1;
    allout1(:,:,cnt1) = out.dfa;
    allpow1(:,:,cnt1) = out.pow;
  else
    cnt2 = cnt2 + 1;
    allout2(:,:,cnt2) = out.dfa;
    allpow2(:,:,cnt2) = out.pow;
  end
  
end


clear out
out.dfa = outa.dfa;
out.pow = outa.pow;

out.dfa_res(:,:,:,1) = allout1(:,:,1:62);
out.dfa_res(:,:,:,2) = allout2(:,:,1:62);

out.pow_res(:,:,:,1) = allpow1(:,:,1:62);
out.pow_res(:,:,:,2) = allpow2(:,:,1:62);

