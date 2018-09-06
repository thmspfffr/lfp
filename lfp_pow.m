% Files heissen mouse_id_s%02d_date_%02d.lfp
% das erste s%02d ist der counter f?r die series / session / penetration
% das zweite %02d ist der counter f?r das experiment (verschiedene Stimuli innerhalb der Session, hier nur Spontanaktivit?t)

fs = 1250;
run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
d = dir('~/lfp/dat/');

freqs = [1:99; 3:101]';
siginfo = nbt_Info;
siginfo.converted_sample_frequency = 1250;
v = 1;

load ~/lfp/expInfo.mat

MOUSE_IDs = {m.mouse_id};
MOUSE_COUNTER = {m.mouse_counter};

outdir = '~/lfp/proc/'

for imouse = 1 : length(MOUSE_IDs)
  
  d = dir(sprintf('~/lfp/dat/%s*.lfp',MOUSE_IDs{imouse}));
  for idir = 1 : length(d)
    
    if ~exist(sprintf([outdir 'lfp_dfa_m%d_d%d_v%d_processing.txt'],imouse,idir,v))
      system(['touch ' outdir sprintf('lfp_dfa_m%d_d%d_v%d_processing.txt',imouse,idir,v)]);
    else
      continue
    end
    
    clear out
    
    out.name = d(idir).name;
    
    for ifoi = 1 : length(freqs)
      
      out.freq = freqs(ifoi,:);

      fid = fopen(['~/lfp/dat/' d(idir).name]);

      lfpData = fread(fid,'int16');
      fclose(fid);

      channel_count = 32;
      % reshape LFP data into channels x samples
      lfpData = reshape(lfpData, channel_count, []);

      ampenv = single(nbt_filter_fir(lfpData',freqs(ifoi,1),freqs(ifoi,2),siginfo.converted_sample_frequency,2/freqs(ifoi,1)));       
      ampenv = abs(hilbert(ampenv));

      tmp = tp_dfa(ampenv,[3 50],fs,0.5,15);
      out.dfa(:,ifoi) = tmp.exp;

    end
    
    save(sprintf([outdir 'lfp_dfa_m%d_d%d_v%d.mat'],imouse,idir,v),'out');
    
  end
end

