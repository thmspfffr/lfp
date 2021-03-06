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
v = 1;
NPERM = 2000;
% --------------------------


mask = find( triu(ones(32)) - eye(32) );

load ~/lfp/expInfo.mat
outdir = '~/lfp/proc/';
% run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m

fs = 1250;

cnt1 = 0; cnt2 = 0;
for idir = 1 : 124
  
  load(['~/lfp/proc/' sprintf('lfp_powcorr_d%d_v%d.mat',idir,v)])
  
  cond(idir)=isempty(strfind(out.mouse_notes,'+/-'));
% 
  if cond(idir) == 1
    cnt1 = cnt1 + 1;
    allout(:,:,:,1,cnt1) = out.powcorr;
    allcoh(:,:,:,1,cnt1) = out.coh;
  else
    cnt2 = cnt2 + 1;
    allout(:,:,:,2,cnt2) = out.powcorr;
    allcoh(:,:,:,2,cnt2) = out.coh;
  end
end

if cnt1~=62 | cnt2~=62
  error('Not all output files500 were loaded!')
end
clear out

for ifoi = 1 : size(allout,3)
  
  [h,p,~,s] = ttest(allout(:,:,ifoi,2,:),allout(:,:,ifoi,1,:),'dim',5,'alpha',0.025);
  out.emp_pow_pos(ifoi) = sum(h(mask) & s.tstat(mask)>0)./length(mask);
  out.emp_pow_neg(ifoi) = sum(h(mask) & s.tstat(mask)<0)./length(mask);

  [h,p,~,s] = ttest(allcoh(:,:,ifoi,2,:),allcoh(:,:,ifoi,1,:),'dim',5,'alpha',0.025);
  out.emp_coh_pos(ifoi) = sum(h(mask) & s.tstat(mask)>0)./length(mask);
  out.emp_coh_neg(ifoi) = sum(h(mask) & s.tstat(mask)<0)./length(mask);
  
end

all_idx = randi(2,[size(allout,5),NPERM]);

for iperm = 1 : NPERM
  iperm
  idx1 = all_idx(:,iperm);
  idx2 = 3-idx1;
  
  for i = 1 : length(idx1)
    
    permdat_pow(:,:,:,1,i) = allout(:,:,:,idx1(i),i);
    permdat_pow(:,:,:,2,i) = allout(:,:,:,idx2(i),i);
    
    permdat_coh(:,:,:,1,i) = allcoh(:,:,:,idx1(i),i);
    permdat_coh(:,:,:,2,i) = allcoh(:,:,:,idx2(i),i);
    
  end
  
  for ifoi = 1 : size(permdat_pow,3)
    
    
    
    [h,p,~,s] = ttest(permdat_coh(:,:,ifoi,2,:),permdat_coh(:,:,ifoi,1,:),'dim',5,'alpha',0.025);
    
    out.perm_coh_neg(ifoi,iperm) = sum(h(mask) & s.tstat(mask)<0)./length(mask);
    out.perm_coh_pos(ifoi,iperm) = sum(h(mask) & s.tstat(mask)>0)./length(mask);
    
    [h,p,~,s] = ttest(permdat_pow(:,:,ifoi,2,:),permdat_pow(:,:,ifoi,1,:),'dim',5,'alpha',0.025);
    
    out.perm_pow_neg(ifoi,iperm) = sum(h(mask) & s.tstat(mask)<0)./length(mask);
    out.perm_pow_pos(ifoi,iperm) = sum(h(mask) & s.tstat(mask)>0)./length(mask);


  end
  
end

out.NPERM  = NPERM;

save(sprintf('~/lfp/proc/lfp_powcorr_permtest_v%d.mat',v),'out');

error('!')

%%
v = 1;

load(['~/lfp/proc/' sprintf('lfp_powcorr_d%d_v%d.mat',1,v)])
freqs = out.freqs;
load(sprintf('~/lfp/proc/lfp_powcorr_permtest_v%d.mat',v));

clear p 

max_emp_pow_pos = max(out.perm_pow_pos);
max_emp_pow_neg = max(out.perm_pow_neg);
max_emp_coh_pos = max(out.perm_coh_pos);
max_emp_coh_neg = max(out.perm_coh_neg);

for ifoi = 1 : size(out.perm_pow_pos,1)
  
  out.p_pow_pos(ifoi) = 1-(sum(out.emp_pow_pos(ifoi) > max_emp_pow_pos)./out.NPERM);
  out.p_pow_neg(ifoi) = 1-(sum(out.emp_pow_neg(ifoi) > max_emp_pow_neg)./out.NPERM);
  out.p_coh_pos(ifoi) = 1-(sum(out.emp_coh_pos(ifoi) > max_emp_coh_pos)./out.NPERM);
  out.p_coh_neg(ifoi) = 1-(sum(out.emp_coh_neg(ifoi) > max_emp_coh_neg)./out.NPERM);

end

save(sprintf('~/lfp/proc/lfp_powcorr_permtest_pval_v%d.mat',v),'out')

%%

figure; hold on
set(gcf,'color','w')

subplot(2,2,1); hold on
plot(out.p_pow_pos,'linewidth',2,'color',[0.9 0.2 0.1])
plot(out.p_pow_neg,'linewidth',2,'color',[0.1 0.2 0.9])
line([find(out.p_pow_pos<0.05,1,'first') find(out.p_pow_pos<0.05,1,'last')],[-0.05 -0.05],'linewidth',3,'color','r')
% line([find(out.p_pow_neg<0.05,1,'first') find(out.p_pow_neg<0.05,1,'last')],[-0.05 -0.05],'linewidth',3,'color','b')
line([1 64],[0.05 0.05],'linestyle',':','color','k')
axis square; axis([1 64 -0.1 1.1]); tp_editplots
set(gca,'XTick',1:9:size(freqs,1),'XTickLabels',mean(freqs(1:9:size(freqs,1),:),2))
xlabel('Frequency [Hz]'); ylabel('P-Value (corrected)')
title('Power correlations')

subplot(2,2,2); hold on
plot(out.p_coh_pos,'linewidth',2,'color',[0.9 0.2 0.1])
plot(out.p_coh_neg,'linewidth',2,'color',[0.1 0.2 0.9])
line([find(out.p_coh_pos<0.05,1,'first') find(out.p_coh_pos<0.05,1,'last')],[-0.05 -0.05],'linewidth',3,'color','r')
line([find(out.p_coh_neg<0.05,1,'first') find(out.p_coh_neg<0.05,1,'last')],[-0.05 -0.05],'linewidth',3,'color','b')
line([1 64],[0.05 0.05],'linestyle',':','color','k')
axis square; axis([1 64 -0.1 1.1]); tp_editplots
set(gca,'XTick',1:9:size(freqs,1),'XTickLabels',mean(freqs(1:9:size(freqs,1),:),2))
xlabel('Frequency [Hz]'); ylabel('P-Value (corrected)')
title('Imaginary coherence')

subplot(2,2,3); hold on
plot(out.emp_pow_pos,'linewidth',2,'color',[0.9 0.2 0.1])
line([find(out.p_pow_pos<0.05,1,'first') find(out.p_pow_pos<0.05,1,'last')],[-0.05 -0.05],'linewidth',3,'color','r')
plot(out.emp_pow_neg,'linewidth',2,'color',[0.1 0.2 0.9])
% line([find(out.p_pow_neg<0.05,1,'first') find(out.p_pow_neg<0.05,1,'last')],[-0.05 -0.05],'linewidth',3,'color','b')
axis square; axis([1 64 -0.1 0.82]); tp_editplots
set(gca,'XTick',1:9:size(freqs,1),'XTickLabels',mean(freqs(1:9:size(freqs,1),:),2))
xlabel('Frequency [Hz]'); ylabel(sprintf('Fraction of \naltered conn. [%%]'))
line([1 64],[0 0 ],'linewidth',1,'color','k')

subplot(2,2,4); hold on
plot(out.emp_coh_pos,'linewidth',2,'color',[0.9 0.2 0.1])
line([find(out.p_coh_pos<0.05,1,'first') find(out.p_coh_pos<0.05,1,'last')],[-0.05 -0.05],'linewidth',3,'color','r')
plot(out.emp_coh_neg,'linewidth',2,'color',[0.1 0.2 0.9])
line([find(out.p_coh_neg<0.05,1,'first') find(out.p_coh_neg<0.05,1,'last')],[-0.05 -0.05],'linewidth',3,'color','b')
axis square; axis([1 64 -0.1 0.82]); tp_editplots
set(gca,'XTick',1:9:size(freqs,1),'XTickLabels',mean(freqs(1:9:size(freqs,1),:),2))
xlabel('Frequency [Hz]'); ylabel(sprintf('Fraction of \naltered conn. [%%]'))
line([1 64],[0 0 ],'linewidth',1,'color','k')

print(gcf,'-dpdf',sprintf('~/lfp/plots/lfp_powcorr_plot_v%d.pdf',v))

