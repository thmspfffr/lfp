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
foi_range       = unique(round(2.^[1:.1:8]));
freqs = [foi_range-(foi_range./2)/2; foi_range+(foi_range./2)/2]';
% --------------------------

load ~/lfp/expInfo.mat
outdir = '~/lfp/proc/';
run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m

fs = 1250;
siginfo = nbt_Info;
siginfo.converted_sample_frequency = fs;

for idir = 1 : length(l)
%   
  if ~exist(sprintf([outdir 'lfp_powcorr_d%d_v%d_processing.txt'],idir,v))
    system(['touch ' outdir sprintf('lfp_powcorr_d%d_v%d_processing.txt',idir,v)]);
  else
    continue
  end

  clear out
  out.freqs = freqs;

  out.exp_num       = e(idir).exp_num;
  out.mouse_counter = e(idir).mouse_counter;
  out.mouse_id      = m([m.mouse_counter]==out.mouse_counter).mouse_id;
  out.series_num    = e(idir).series_num;
  out.mouse_notes   = m([m.mouse_counter]==out.mouse_counter).mouse_notes;
  out.mean_speed    = l(idir).mean_speed;

  dd = dir(['~/lfp/dat/' sprintf('%s_s%02d_*_%02d.lfp',out.mouse_id, out.series_num,out.exp_num)]);
  
  fid = fopen(['~/lfp/dat/' dd.name]);
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
   
  save(sprintf([outdir 'lfp_powcorr_d%d_v%d.mat'],idir,v),'out');
    
end

error('!')

%% NEED TO COMPUTE PERMTEST FOR RESULT THAT CAN BE PLOTTED IN AN EASY WAY

v = 1;
clear allpow allout allcoh
load ~/lfp/expInfo.mat
addpath('/home/tpfeffer/Documents/MATLAB/Colormaps/Colormaps (5)/Colormaps/')
dd = dir(sprintf('~/lfp/proc/*powcorr*v%d*mat',v));

allpowcorr   = nan(40,40,size(freqs,1),2,length(dd)/2);
allimagcoh   = nan(40,40,size(freqs,1),2,length(dd)/2);
all_depth = [min([d(:).rel_depth]):25:max([d(:).rel_depth])];

cnt1 = 0; cnt2 = 0;
for idir = 1 : length(dd)
  
  load(['~/lfp/proc/' sprintf('lfp_powcorr_d%d_v%d.mat',idir,v)])
  % ALIGN RECORDINGS
  depth_idx       = [d.mouse_counter]==out.mouse_counter&[d.series_num]==out.series_num;
  rel_depth       = min([d(depth_idx).rel_depth]):25:max([d(depth_idx).rel_depth]);
  indiv_depth_idx = ~isnan(d(depth_idx).rel_depth);
  abs_depth_idx   = all_depth>=rel_depth(1) & all_depth<=rel_depth(end);
    
  cond(idir)=isempty(strfind(out.mouse_notes,'+/-'));
% 
  if cond(idir) == 1
    cnt1 = cnt1 + 1;
    allpowcorr(abs_depth_idx,abs_depth_idx,:,1,cnt1) = out.powcorr(indiv_depth_idx,indiv_depth_idx,:);
    allimagcoh(abs_depth_idx,abs_depth_idx,:,1,cnt1) = out.coh(indiv_depth_idx,indiv_depth_idx,:);
  else
    cnt2 = cnt2 + 1;
    allpowcorr(abs_depth_idx,abs_depth_idx,:,2,cnt2) = out.powcorr(indiv_depth_idx,indiv_depth_idx,:);
    allimagcoh(abs_depth_idx,abs_depth_idx,:,2,cnt2) = out.coh(indiv_depth_idx,indiv_depth_idx,:);
  end
  
  freqs = out.freqs;
end

if cnt1~=62 | cnt2~=62
  warning('Not all output files were loaded.')
end


out.dfa = allpowcorr;
out.coh = allcoh;


  
  
  


