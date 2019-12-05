function codes = stimcodes_FST

% This file contains the 'timestamp' and stimulus 'tag' codes that are
% Strobed by PLDAPS to the Omniplex system.
% WARNING: This is the key to the data. DO NOT change the contents of this
% file.

%% 
% trial codes
codes.trialBegin = 30001;%1001;
codes.trialEnd = 30009;%1009;

codes.connectPLX = 11001;
codes.trialcount = 11002;
codes.blocknumber = 11003;
codes.trinblk     = 11004;
codes.setnumber   = 11005;
codes.state       = 11008;
codes.trialcode   = 11009;
codes.trialtype   = 11010;
codes.tasktype    = 11099; 
% 1 = Visually guided saccades; 2 = Memory guided saccades; 3 = Attn Motion
% task; 4 = Microstimulation; 5- dir tuning; 6- Orn tuning; 7 - Attn Orn
% Task; 8 - RF mapping
codes.repeat20 = 11098; % 1 = 20 repeat trials during MemSac task.
codes.vissac    = 11097; % 1 = vis sac; 0 = memsac protocol
codes.inactivation = 11095; % during inactivation

% joystick codes
codes.joypress = 2001;
codes.joyrelease = 2002;
codes.joybreak = 2005;

% fixation codes
codes.fixdoton = 3001;
codes.fixdotdim = 3002;
codes.fixdotoff = 3003;
codes.fixacqd = 3004;
codes.fixbreak = 3005;

codes.fixTH = 13001;
codes.fixR = 13002;
codes.fixdotdimvalue = 13003;
codes.fixchangetrial = 13004;

% saccade codes (used in VisSac and MemSac tasks)
codes.saccadeonset = 2003;
codes.saccadeoffset = 2004;

% target codes (used in VisSac and MemSac tasks)
codes.targeton = 4001;
codes.targetdim = 4002;
codes.targetoff = 4003;
codes.targetacqd = 4004;
codes.targfixbreak = 4005;

codes.targetTH = 14001;
codes.targetR = 14002;

% cue codes (used in Attn tasks)
codes.cueon = 5001;
codes.cueoff = 5003;

codes.cueTH = 15001;
codes.cueR = 15002;
codes.cuecolor = 15003;

% stimulus codes (used in Attn tasks)
codes.stimon        = 6002;
codes.stimoff       = 6003;
codes.cuechange     = 6004;
codes.foilchange    = 6005;
codes.nochange      = 6006;

codes.RFlocecc      = 16001;
codes.RFlocTH       = 16002;

codes.stimchangetrial = 16003;
codes.changeloc     = 16004;

% codes for PA motion task
codes.loc1dir       = 16005;
codes.loc2dir       = 16006;
codes.loc1del       = 16007;
codes.loc2del       = 16008;

% codes for PA orientation task
codes.loc1orn       = 16005;
codes.loc2orn       = 16006;
codes.loc1amp       = 16007;
codes.loc2amp       = 16008;

% code for motion directions in dir tuning task
codes.motdir = 24000; % this is to send stim info; for eevnt codes for each dir see trialcodes.dirtun.m

% code for orientations in orn tuning task
codes.orn = 25000; % this is to send stim info; for eevnt codes for each orn see trialcodes.orntun.m

% % code for objects in obj tuning task
% codes.obj = 26000; % this is to send stim info; for eevnt codes for each obj see trialcodes.objtun.m

% reward code
codes.reward = 8000;
codes.rewardduration = 18000;

% micro stim codes
codes.microstimon = 7001;
codes.saconset = 7005;

% audio codes
codes.audioFBKon = 9000;
codes.lowtone = 9001;
codes.noisetone = 9002;
codes.hightone = 9003;

