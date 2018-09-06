%% lfp_powcorr
% analyses lfps recorded from several layers across mouse V1
% during periods of locomotion and periods of stationarity. 
% Output generated is in 'out' structure. Contains prt, power and
% amplitudes computed across various frequencies of interest.

% Next script: lfp_powcorr_permtest.m

clear

% --------------------------
% VERSION 1
% --------------------------
v = 2;
foi_range       = unique(round(2.^[1:.1:9]));
freqs = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
freqs = freqs(1:end-1,:);
% --------------------------

load ~/lfp/expInfo.mat
outdir = '~/lfp/proc/';
run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m

fs = 1250;
siginfo = nbt_Info;
siginfo.converted_sample_frequency = fs;

for idir = 1 : length(l)
  
  if ~exist(sprintf([outdir 'lfp_stim_powcorr_d%d_v%d_processing.txt'],idir,v))
    system(['touch ' outdir sprintf('lfp_stim_powcorr_d%d_v%d_processing.txt',idir,v)]);
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
  
  d = dir(['~/lfp/dat/' sprintf('%s_s%02d_*_%02d.lfp',out.mouse_id, out.series_num,out.exp_num)]);
  
  fid = fopen(['~/lfp/dat/' d.name]);
  lfpData = fread(fid,'int16');
  fclose(fid);
  
  channel_count = 32;
  % reshape LFP data into channels x samples
  lfpData = reshape(lfpData, channel_count, []);
  cs = data2cs_event(lfpData',2500,1250,size(lfpData,2),[]);
  cs = cs(:,:,2:end);
  f  = 0.5:0.5:fs/2;
  
  for ii = 1 : size(cs,3)
    tmp_coh(:,:,ii) = imag(cs(:,:,ii)./sqrt(diag(cs(:,:,ii))*diag(cs(:,:,ii))'));
  end
  
  for ifoi = 1 : length(freqs)
    fprintf('Computing icoh for freq%d...\n',ifoi)

    idx = f>=freqs(ifoi,1)&f<=freqs(ifoi,2);
    out.coh(:,:,ifoi) = mean(tmp_coh(:,:,idx),3);  
    
    fprintf('Computing powcorr for freq%d...\n',ifoi)
    ampenv = nbt_filter_fir(lfpData',freqs(ifoi,1),freqs(ifoi,2),fs,2/freqs(ifoi,1));
    ampenv = hilbert(ampenv);
    out.powcorr(:,:,ifoi) = compute_orthopowcorr(ampenv');
  
  end
  
  save(sprintf([outdir 'lfp_stim_powcorr_d%d_v%d.mat'],idir,v),'out');
    
end

error('!')

%%
v = 2;
addpath('/home/tpfeffer/Documents/MATLAB/Colormaps/Colormaps (5)/Colormaps/')

% d = dir(sprintf('~/lfp/proc/*powcorr*v%d*mat',v));
clear allpow allout allcoh

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

v = 2;
addpath('/home/tpfeffer/Documents/MATLAB/Colormaps/Colormaps (5)/Colormaps/')
clear allpow allout allcoh

cnt1 = 0; cnt2 = 0;
for idir = 1 : 124
  
  load(['~/lfp/proc/' sprintf('lfp_powcorr_d%d_v%d.mat',idir,v)])
  
  cond(idir)=isempty(strfind(out.mouse_notes,'+/-'));
% 
  if cond(idir) == 1
    cnt1 = cnt1 + 1;
    allout(:,:,:,1,cnt1) = out.powcorr;
    allcoh(:,:,:,1,cnt1) = out.coh;
%     allpow(:,:,1,cnt1) = out.pow;
%     alldepth(cnt1,1) = depth;
  else
    cnt2 = cnt2 + 1;
    allout(:,:,:,2,cnt2) = out.powcorr;
    allcoh(:,:,:,2,cnt2) = out.coh;
%     allpow(:,:,2,cnt2) = out.pow;
%     alldepth(cnt2,2) = depth;
  end
  
  freqs = out.freqs;
end

if cnt1~=62 | cnt2~=62
  warning('Not all output files were loaded.')
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

% close;

if exist(sprintf('~/lfp/proc/lfp_powcorr_permtest_pval_v%d.mat',v))
  load(sprintf('~/lfp/proc/lfp_powcorr_permtest_pval_v%d.mat',v))
  outp = out;
end

out.dfa = allout;
out.coh = allcoh;

%%
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
  imagesc(par,[-0.15 0.15]); colormap(parula)
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
tp_editplots
set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',[])

subplot(2,2,3)

subplot(4,4,[9 10 13 14])
[~,p] = ttest(out.dfa(:,:,2,:),out.dfa(:,:,1,:),'dim',4);
% p_mask]=fdr(p,0.2);
par = mean(out.dfa(:,:,2,:),4)-mean(out.dfa(:,:,1,:),4);
% par(p_mask==0)=nan;
imagescnan(-log10(p),[0 5]); colormap(parula)
xlabel('Frequency [Hz]'); ylabel('Channel #')

set(gca,'YTick',1:5:32,'YTickLabels',1:5:32)
set(gca,'XTick',1:9:size(out.freqs,1)-1,'XTickLabels',mean(out.freqs(1:9:end-1,:),2))
tp_editplots

subplot(4,4,[11 15])
par = flipud(mean(mean(out.dfa(:,:,2,:),4),2))-flipud(mean(mean(out.dfa(:,:,1,:),4),2));
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
%%

out.freqs = freqs;
ifoi = out.p_pow_pos<0.05
close
figure; set(gcf,'color','w')
subplot(2,2,1);
par = abs(out.dfa(:,:,ifoi,:,:));
[~,p,~,s]=ttest(mean(par(:,:,:,2,:),3),mean(par(:,:,:,1,:),3),'dim',5);
p_corr = fdr(p); 
p_mask = (p_corr) < 0.05 & s.tstat>0;
t = s.tstat;
t(p_mask==0)= nan; 
imagescnan(t,[-6 6]); axis square
colormap(parula)
xlabel('Channel #'); ylabel('Channel #')
title('oPEC: Pos')
tp_editplots

subplot(2,2,2);
ifoi = out.p_pow_neg<0.05
par = abs(out.dfa(:,:,ifoi,:,:));
[~,p,~,s]=ttest(mean(par(:,:,:,2,:),3),mean(par(:,:,:,1,:),3),'dim',5);
p_corr = fdr(p); 
p_mask = (p_corr) < 0.05 & s.tstat<0;
t = s.tstat;
t(p_mask==0)= nan; 
imagescnan(t,[-6 6]); axis square
colormap(parula)
xlabel('Channel #'); ylabel('Channel #')
title('oPEC: Neg')
tp_editplots

ifoi = out.p_coh_pos<0.05
subplot(2,2,3);
par = abs(out.coh(:,:,ifoi,:,:));
[~,p,~,s]=ttest(mean(par(:,:,:,2,:),3),mean(par(:,:,:,1,:),3),'dim',5);
p_corr = fdr(p); 
p_mask = (p_corr) < 0.05 & s.tstat>0;
t = s.tstat;
t(p_mask==0)= nan;

imagescnan(t,[-6 6]); axis square
colormap(parula)

xlabel('Channel #'); ylabel('Channel #')
title('iCoherence: Pos.')

ifoi = out.p_coh_neg<0.05
subplot(2,2,4);
par = abs(out.coh(:,:,ifoi,:,:));
[~,p,~,s]=ttest(mean(par(:,:,:,2,:),3),mean(par(:,:,:,1,:),3),'dim',5);
p_corr = fdr(p); 
p_mask = (p_corr) < 0.05 & s.tstat<0;
t = s.tstat;
t(p_mask==0)= nan;

imagescnan(t,[-6 6]); axis square
colormap(parula)
title('iCoherence: Neg.')
xlabel('Channel #'); ylabel('Channel #')

print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_powcorr_contrast_v%d.pdf',v))

%%
figure;
subplot(2,2,3); hold on
for ifoi = 1 : size(out.freqs,1)
  
  [~,p,~,s]=ttest(abs(out.dfa(:,:,ifoi,2,:)),abs(out.dfa(:,:,ifoi,1,:)),'dim',5);
  p_corr = fdr(p); mask = p_corr<0.05; 
  
  ss1(ifoi) = sum(s.tstat(mask)>0)./length(p_corr(:));
  ss2(ifoi) = sum(s.tstat(mask)<0)./length(p_corr(:));
end

plot(ss1,'linewidth',3,'color',[1 0.5 0.2]); hold on
plot(ss2,'linewidth',3,'color',[0.2 0.5 1])

axis square
xlabel('Frequency [Hz]'); ylabel('Fraction of alt. conn.')
axis([0 64 -0.1 0.6])
set(gca,'xtick',1:8:64,'xticklabel',mean(out.freqs(1:8:end,:),2))

subplot(2,2,4); hold on
for ifoi = 1 : size(out.freqs,1)
  
  [~,p,~,s]=ttest(abs(out.coh(:,:,ifoi,2,:)),abs(out.coh(:,:,ifoi,1,:)),'dim',5);
  p_corr = fdr(p); 
  ss1(ifoi) = sum(s.tstat(mask)>0)./length(p_corr(:));
  ss2(ifoi) = sum(s.tstat(mask)<0)./length(p_corr(:));
  
end

plot(ss1,'linewidth',3,'color',[1 0.5 0.2]); hold on
plot(ss2,'linewidth',3,'color',[0.2 0.5 1])

axis square
xlabel('Frequency [Hz]'); ylabel('Fraction of alt. conn.')
axis([0 64 -0.1 0.6])
set(gca,'xtick',1:8:64,'xticklabel',mean(out.freqs(1:8:end,:),2))

%% CORRELATION POWER CORRELATION / IMAGINARY COHERENCE
mask = find( triu(ones(32)) - eye(32) );

for ifoi = 1 : size(out.freqs,1)
  c1 = squeeze(out.coh(:,:,ifoi,1,:));
  c2 = squeeze(out.dfa(:,:,ifoi,1,:));
  
  for isubj = 1 : 60
    tmp1 = c1(:,:,isubj); tmp1 = tmp1(mask);
    tmp2 = c2(:,:,isubj); tmp2 = tmp2(mask);
    
    c(ifoi,isubj) = corr(tmp1(:),tmp2(:));
  end
end

    
    
  
  
  
  


