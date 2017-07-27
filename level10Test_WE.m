function level10Test_WE

% White Elephant Level 8:
% Two or more stimuli are presented
% Presentations are limited to hold times (i.e. no repetition)
% Mistakes lead to time outs

global DA gf h saveData
% DA: TDT connection structure
% gf: Go ferrit user data
% h:  Online GUI handles

%GUI Clock
gf.sessionTime = (now - gf.startTime)*(24*60*60);
set(h.sessionLength,'string', sprintf('%0.1f secs',gf.sessionTime));

gf.level = 10;

% Reverse gf.rewardTimes from parameters
% It's [centre, 2, 3, 4, 5, 6, 7, 8]
% Should be [8, 7, 6, 5, 4, 3, 2, centre]
gf.rewardTimes = gf.rewardTimes(end:-1:1);

%Run case

switch gf.status
    
    % __________________________________________________________________________
    case('PrepareStim')
        
        gf.response=[];
        
        % Identify as new trial
        gf.correctionTrial = 0;
        
        calculateHoldTimes
        calculateAtten
        
        % Create 12 row zeros matrix for sounds here
        stimSound = zeros(12,ceil(gf.cutOff/1000*gf.fStim));
        stimLight = zeros(12,ceil(gf.cutOff/1000*gf.fStim));
        
        % Select trial type
        gf.TT=gf.trialTypes(randi(numel(gf.trialTypes)));
        % stimSound=rand(gf.nRefSpks+1,floor(gf.fStim)*3);
        % Create 7 row zeros matrix for sounds here
        stimSound = zeros(12,ceil(gf.cutOff/1000*gf.fStim));
        soundLR = zeros(ceil(gf.cutOff/1000*gf.fStim),2);
        stimLight = zeros(12,ceil(gf.cutOff/1000*gf.fStim));
        lightLR = zeros(ceil(gf.cutOff/1000*gf.fStim),2);
        
        % soundLR and light LR will contain L and R sounds
        % These will be written to stimSound and stimLight (row for each
        % chan)
        % Then whole of stimSound and stimLight written to buffer to make
        % correct channels are blanked

        % Generate stims for active channels here
        
        % Select combination of events
        gf.evInd = randi(numel(gf.eventsFast));
        
        events = round(rand)+1; % should really be called side
        if events==1
            eV = [gf.eventsSlow(gf.evInd), gf.eventsFast(gf.evInd)];
            gf.side=2;
            disp('Answer should be left')
        else
            eV = [gf.eventsFast(gf.evInd), gf.eventsSlow(gf.evInd)];
            gf.side=8;
            disp('Answer should be right')
        end
        gf.tarIdx=[2,8];
        disp(eV(end:-1:1)) % NB eV 1=right 2=left!
        
        % Generates soundLR(:,[L,R]) and lightLR(:,[L,R])
        % (event for each determined above)
        % Pick noise to use on all channels
        % (Picks for both modalities, regardless of type)
        gf.nA = gf.aNoise(randi(numel(gf.aNoise)));
        cfg.noiseMag=gf.nA;
        gf.nV = gf.vNoise(randi(numel(gf.vNoise)));
        cfg.noiseMag=gf.nV;
        for i=1:2 % Generating for 2 sides
            clear cfg
            
            % Common parameters
            % p=gf.speakers(i);
            % (here left=8 and right=2, so ev processed in correct order)
            cfg.nEvents = eV(i);
            cfg.cutOff = gf.cutOff;
            cfg.eventMag = gf.eventMag;
            cfg.gap1 = gf.gap1;
            cfg.gap2 = gf.gap2;
            cfg.Fs = gf.fStim;
            cfg.eventType = gf.eventType;
            cfg.rise = gf.stimRise;
            cfg.duration = gf.duration;
            cfg.noiseType = 'multipleBlocks';
            cfg.rideNoise = gf.rideNoise;
            cfg.dispWarn = 0;
            
            % Specific parameters and generation
            gf.aysncOffset=[];
            switch gf.TT
                case 2 % AO
                    cfg.type='Aud';
                    cfg.rise=gf.stimRise;
                    % cfg.noiseMag=gf.aNoise;
                    cfg.engBuff=gf.endBuff;
                    cfg.startBuff=gf.startBuff;
                    cfg.rideNoise=gf.rideNoise;
                    
                    % Get noise
                    cfg.noiseMag=gf.nA;

                    % Seed
                    cfg.seed=gf.seeds(1,1);
                    
                    aStim=Church2(cfg);
                    vStim=NaN;
                    
                case 3 % VO
                    cfg.eventType='flat';
                    cfg.eventMag=gf.vRange;
                    cfg.engBuff=gf.endBuff;
                    cfg.startBuff=gf.startBuff;
                    % cfg.noiseMag=gf.vNoise;
                    cfg.rideNoise=gf.rideNoise;
                    % Get noise
                    cfg.noiseMag=gf.nV;
                    % Seed
                    cfg.seed=gf.seeds(2,2);
                    
                    aStim=NaN;
                    vStim=Church2(cfg);
                    
                case 4 % AV sync
                    cfg.type='Aud';
                    cfg.rise=gf.stimRise;
                    cfg.engBuff=gf.endBuff;
                    cfg.startBuff=gf.startBuff;
                    cfg.noiseMag=gf.aNoise;
                    cfg.rideNoise=gf.rideNoise;
                    % Seed
                    cfg.seed=gf.seeds(3,1);
                    
                    aStim=Church2(cfg);
                    
                    % Keep cfg, replace V fields
                    cfg.type='Vis';
                    cfg.eventType='flat';
                    cfg.noiseMag=gf.vNoise;
                    cfg.eventMag=gf.vRange;
                    % Seed
                    cfg.seed=gf.seeds(3,2); % (Same as A)
                    
                    vStim=Church2(cfg);
                    
                case 5 % AV async
                    cfg.type='Aud';
                    cfg.engBuff=gf.endBuff;
                    cfg.startBuff=gf.startBuff;
                    cfg.noiseMag=gf.aNoise;
                    cfg.rideNoise=gf.rideNoise;
                    % Seed
                    cfg.seed=gf.seeds(4,1);
                    
                    aStim=Church2(cfg);
                    
                    clear cfg
                    
                    % Re-set common
                    cfg.Fs = gf.fStim; % Sampling rate
                    cfg.eventLength = gf.stimEventDuration;
                    %                     cfg.mag = 0-gf.atten;
                    cfg.eventType = gf.eventType;
                    cfg.gap1 = gf.gap1; % Gap durations
                    cfg.gap2 = gf.gap2;
                    cfg.duration=gf.duration;
                    cfg.cutOff=gf.cutOff;
                    cfg.eventFreq = gf.eventFreq; % If 'sine', frequency
                    cfg.noiseType='multipleBlocks';
                    cfg.rise=gf.stimRise;
                    cfg.rideNoise=gf.rideNoise;
                    
                    % Select variable offset
                    gf.aysncOffset = ...
                        gf.asyncOffsets(randi(numel(gf.asyncOffsets)));
                    
                    % Set specific
                    cfg.type='Vis';
                    cfg.eventType='flat';
                    cfg.eventMag=gf.vRange;
                    cfg.endBuff=gf.endBuff+gf.aysncOffset;
                    cfg.startBuff=gf.startBuff+gf.aysncOffset;
                    cfg.noiseMag=gf.vNoise;
                    % Seed
                    cfg.seed=gf.seeds(4,2);
                    
                    vStim=Church2(cfg);
            end
            
            % Stimrecord is : rows = mod, cols = side (1=L, 2=R)
            if isstruct(aStim)
                soundLR(1:length(aStim.sound),i) = ...
                    aStim.sound';
                aStim = rmfield(aStim, 'sound');
                gf.stimRecord{1,i}=aStim;
                gf.stimRecord{2,i}=NaN;
            end
            if isstruct(vStim)
                lightLR(1:length(vStim.sound),i) = ...
                    (vStim.sound*gf.vMulti)';
                vStim = rmfield(vStim, 'sound');
                gf.stimRecord{2,i}=vStim;
                gf.stimRecord{1,i}=NaN;
            end
        end

        % 22/2/2016
        % Select random row from gf.speakers
        gf.sr = randi(size(gf.speakers,1));
        % Copy sounds for each side to relevant channels
        % spks 2-4 = R, 5= centre, 6-8=L
        for i = 1:length(gf.speakers)
            p = gf.speakers(gf.sr,:);
            if p(i)==5
                % Centre position
            elseif any(p(i)==[6,7,8]) % Left positions
                stimSound(p(i),:) = soundLR(:,1)';
                stimLight(p(i),:) = lightLR(:,1)';
            elseif any(p(i)==[2,3,4]) % Right positions
                stimSound(p(i),:) = soundLR(:,2);
                stimLight(p(i),:) = lightLR(:,2)';
            end
        end
        
        % Apply the attenuation
        stimSound=stimSound.* 10^(-(gf.atten/20));
        
        % Write sound to buffers
        for i=2:8
            DA.WriteTargetVEX(['RX8.speaker', num2str(i)], ... 
                0, 'F32', stimSound(i,:));
        end
        for i=2:8
            DA.WriteTargetVEX(['RX8.ledin', num2str(i)], 0, ...
                'F32', stimLight(i,:)) % it's really light in level 5
        end
        
        % Make holdOK length of entire sound
        % Total hold time is hold time + holdOK
        holdOK = ... 
            ((gf.refDuration*gf.nRefSpks)+(gf.tarDuration/2))*gf.fStim;                % point in stimulus that the ferret must hold to
        % holdOK      = length(stimSound)
        %         for ii=2:8
        %             silentSpks=strcat('%s.speaker',num2str(ii));
        %             silencefiller=DA.WriteTargetVEX(sprintf(silentSpks, gf.stimDevice), 0, 'F32', silence); % fill unused speakers with silence
        %             %  end
        %         end
        
        % Calculate timing information
        playDelay = gf.holdSamples - holdOK;
        %playDelay =0;
        refractS = playDelay +length(stimSound) ...
            + ceil(gf.refractoryTime * gf.fStim);
        absentS = ceil(gf.absentTime * gf.fStim);
        
        % Enable spout LEDs
        DA.SetTargetVal( sprintf('%s.LEDthreshold', ...
            gf.stimDevice), 0.01);
        
        % Set length of time to play noise during timeout
        DA.SetTargetVal( sprintf('%s.noiseDuration', ...
            gf.stimDevice), gf.noiseDuration);
        
        % Enable / Disable Circuit components
        DA.SetTargetVal( sprintf('%s.repeatPlayEnable', gf.stimDevice), ...
            0);                 % Disable OpenEx driven sound repetition
        DA.SetTargetVal( sprintf('%s.repeatPlay',       gf.stimDevice), ...
            0);                 % Disable Matlab driven sound repetition
        DA.SetTargetVal( sprintf('%s.repeatPlayEnable', gf.stimDevice), ...
            0);                 % Disable OpenEx driven sound repetition
        
        % Set timing information on TDT
        DA.SetTargetVal( sprintf('%s.stimNPoints', gf.stimDevice), ...
            size(stimSound,2));
        DA.SetTargetVal( sprintf('%s.firstNPoints', gf.stimDevice), ...
            (size(stimSound,2) + gf.repInterval*gf.fStim)/gf.fStim*1000);
        DA.SetTargetVal( sprintf('%s.repNPoints', gf.stimDevice), ...
            size(stimSound,2) + gf.repInterval*gf.fStim);
        DA.SetTargetVal( sprintf('%s.holdSamples', gf.stimDevice), ...
            gf.holdSamples);
        DA.SetTargetVal( sprintf('%s.absentSamps', gf.stimDevice), ...
            absentS);
        DA.SetTargetVal( sprintf('%s.playDelay', gf.stimDevice), ...
            playDelay);
        DA.SetTargetVal( sprintf('%s.refractorySamps', gf.stimDevice), ...
            refractS);
        
        
        % Enable / Disable Circuit components
        DA.SetTargetVal( sprintf('%s.centerEnable',  ...
            gf.stimDevice), 1);
        if gf.repeatMode==0
            DA.SetTargetVal(sprintf('%s.repeatPlayEnable', ...
                gf.stimDevice), 0); % Disable OpenEx driven sound repetition
            DA.SetTargetVal(sprintf('%s.repeatPlay', ...
                gf.stimDevice), 0); % Disable Matlab driven sound repetition
            DA.SetTargetVal(sprintf('%s.repeatPlayEnable', ...
                gf.stimDevice), 0); % Disable OpenEx driven sound repetition
            DA.SetTargetVal(sprintf('%s.RPMon', gf.stimDevice), 0);
        elseif gf.repeatMode==1
            DA.SetTargetVal(sprintf('%s.repeatPlayEnable', gf.stimDevice), 0);
            % Doesn't need to be on to repeat?
            DA.SetTargetVal(sprintf('%s.repeatPlay', gf.stimDevice), 0);
            % RPMon turned on when trial activated, to avoid lights being
            % left on with centre lick but fialed trigger
            DA.SetTargetVal(sprintf('%s.RPMon', gf.stimDevice), 0);
        end
        
        % Update online GUI
        set(h.status,     'string',sprintf('%s',gf.status))
        set(h.side,       'string',num2str( gf.side))
        set(h.pitch,      'string','N/A')
        set(h.holdTime,   'string',sprintf('%.0f ms', gf.holdTime))
        % set(h.currentStim,'string',num2str(gf.tarChoice))
        set(h.atten,      'string',sprintf('%.1f dB', gf.atten))
        set(h.trialInfo,  'string',sprintf('%d', gf.TrialNumber - 1))
        
        
        gf.status='WaitForStart';
        
        % Center Response__________________________________________________
    case('WaitForStart')
        
        DA.SetTargetVal( sprintf('%s.LEDthreshold',     gf.stimDevice), 0.01);
        DA.SetTargetVal( sprintf('%s.RPMon',       gf.stimDevice), 0);

        % Reset the timeout player
        gf.timeoutPlayed=0;
       
        
        DA.SetTargetVal( sprintf('%s.ledEnable', gf.stimDevice), 1);                 % Enable constant LED in hold time
        DA.SetTargetVal( sprintf('%s.spoutPlayEnable', gf.stimDevice), 1);                 % Enable sound in hold time
        %
        centerLick  = invoke(DA,'GetTargetVal',sprintf('%s.lick1',gf.stimDevice));

        
        if centerLick == 0; % If no start
            % Flash LED
            DA.SetTargetVal(sprintf('%s.led1in',gf.stimDevice),1);
            comment = 'LED flashing, waiting for center lick';
        else % Start presentation
            DA.SetTargetVal( sprintf('%s.led1in', gf.stimDevice), 0);
            DA.SetTargetVal( sprintf('%s.repeatPlayEnable', ...
                gf.stimDevice), 0); % Disable OpenEx driven sound repetition
            
            if gf.repeatMode==1
                DA.SetTargetVal( sprintf('%s.repeatPlay', ...
                    gf.stimDevice), 1); % Disable Matlab driven sound repetition
                DA.SetTargetVal( sprintf('%s.RPMon', ...
                    gf.stimDevice), 1) % Repeat triggered, LEDEnable2 disabled
            end

            gf.startTrialTime = ...
                DA.GetTargetVal(sprintf('%s.lick1Time',gf.stimDevice)) ...
                ./ gf.fStim;  %Open Ex
            gf.status = 'WaitForResponse';
            
            % Reward at center spout
            if gf.centerRewardP > rand(1),
                gf.centerReward = 1;
                comment         = 'Center spout licked - giving reward';
                valve_WE(1,100,1)
            else
                gf.centerReward = 0;
                comment         = 'Center spout licked - no reward';
            end
        end
        
        %Update GUI
        set(h.status,'string',gf.status);
        set(h.comment,'string',comment);
        % disp(['Tick 5:', num2str(toc)])
        %         end
        
        % Peripheral Response______________________________________________________
    case('WaitForResponse')
        tic
        if gf.useLEDs,
            DA.SetTargetVal( sprintf('%s.LEDthreshold',     gf.stimDevice), 0.01);          % Enable LEDs
        end
       
        % Check for peripheral response
        pResp = zeros(2,1);
        gf.tarIdx=[2,8];
        
        for ii = 1 : 2
            pResp(ii) = DA.GetTargetVal( sprintf( '%s.lick%d',...
                gf.stimDevice, gf.tarIdx(ii)));
        end
        % disp(['Tick 2:', num2str(toc)])
        
        % Separate response types
        errResp = pResp;
        % Remove index of correct response
        errResp(gf.tarIdx==gf.side) = [];
        corrResp = pResp(gf.tarIdx==gf.side);
        
        % If no response
        if ~sum(pResp)
            timeNow = ...
                DA.GetTargetVal(sprintf('%s.zTime',gf.stimDevice)) ...
                ./ gf.fStim;
            timeElapsed = timeNow - gf.startTrialTime;
            timeRemaining = gf.abortTrial - timeElapsed;
            
            comment = ...
                sprintf('Awaiting response: \nTime remaining %0.1f s', ...
                timeRemaining); %#ok<NASGU>

            % Check response countdown
            if timeRemaining <= 0 % Abort
                
                % Disable LEDs
                DA.SetTargetVal(sprintf('%s.LEDthreshold', ...
                    gf.stimDevice), 99);
                
                % Log aborted response
                gf.responseTime = -1;
                response = -1;
                
                gf.trialLog(:,1) = [gf.TrialNumber; gf.side; response; corrResp];
                gf.TrialNumber  = gf.TrialNumber + 1;
                % Remove filters - slow to load and save
                saveData{1,gf.TrialNumber} = rmfield(gf, 'SpkrLowPass');
                save([gf.filename2(1:end-4),'.mat'],'saveData');
                
                % Update perfomance graph
                updatePerformance(3) % code 3 = abort trial
                
                gf.status = 'PrepareStim';
            end
            
         else % If response
            % Disable LEDs
            DA.SetTargetVal(sprintf('%s.LEDthreshold', gf.stimDevice), 99);              
                        
            % If animal responds correctly
            if corrResp
                DA.SetTargetVal(sprintf('%s.repeatPlayEnable', ...
                    gf.stimDevice), 0); % Disable OpenEx driven sound repetition
                DA.SetTargetVal(sprintf('%s.repeatPlay', ...
                    gf.stimDevice), 0);
                
                if gf.side == 2;
                    gf.resp = 2;
                elseif gf.side == 8
                    gf.resp = 8;
                end
                
                %Log response
                gf.responseTime = ...
                    DA.GetTargetVal(sprintf('%s.lick%dtime', ...
                    gf.stimDevice, gf.side)) ./ gf.fStim;
                
                gf.LeftRight=[1,0]; % reversed because if tarChoice =2 we need response to the right (1) not left (0)
                response = gf.LeftRight(gf.tarIdx==gf.side);

                % Log
                gf.trialLog(:,1) = ...
                    [gf.TrialNumber; gf.side; response; corrResp];
                
                gf.TrialNumber  = gf.TrialNumber + 1;
                saveData{1,gf.TrialNumber} = rmfield(gf, 'SpkrLowPass');
                save([gf.filename2(1:end-4),'.mat'],'saveData');
                
                % 17/10/2016
                % Reward based on specific spout time
                valve_WE(gf.resp, gf.rewardTimes(gf.resp), 1)
                comment = sprintf('Correct response at speaker %d - giving reward', gf.side);
                
                gf.status  = 'PrepareStim';
                
                % Update perfomance graph
                if response==1
                    updatePerformance(2)             % code 2 = right correct
                elseif response==0
                    updatePerformance(4)             % code 4 = left correct
                end
                
                % If animal responds incorrectly
            elseif any(errResp)
                gf.timeoutPlayed=0;
                gf.LeftRight=[1,0]; % reversed because if tarChoice =2 we need response to the right (1) not left (0)
                if gf.side==2;
                    gf.resp=8;
                elseif gf.side==8
                    gf.resp=2;
                end
                response        = gf.LeftRight(gf.tarIdx~=gf.side);
                
                %%%% ST Modification (Added to allow code to run: 26 Nov
                %%%% 2012 AM.
                %Log response
                gf.responseTime = DA.GetTargetVal(sprintf('%s.lick%dtime',gf.stimDevice, gf.tarIdx(gf.tarIdx~=gf.side))) ./ gf.fStim;
                
                comment    = sprintf( 'Incorrect response at speaker %d', gf.tarIdx(pResp == 1));
                
                
                % logTrial(gf.centerReward, response)                   %See toolbox
                gf.trialLog(:,1) = [gf.TrialNumber; gf.side; response; corrResp];
                gf.TrialNumber  = gf.TrialNumber + 1;
                saveData{1,gf.TrialNumber} = rmfield(gf, 'SpkrLowPass');
                save([gf.filename2(1:end-4),'.mat'],'saveData');
                
                % Update perfomance graph
                if response==1
                    updatePerformance(1)             % code 1 = right incorrect
                elseif response==0
                    updatePerformance(5)             % code 5 = left incorrect
                end
                
                
                % Disable center spout
                DA.SetTargetVal( sprintf('%s.centerEnable', gf.stimDevice), 0);
                
                gf.status  = 'timeout';
                
            end
            
            
        end
        
        % Timeout _________________________________________________________________
        
    case('timeout')
        
        DA.SetTargetVal( sprintf('%s.repeatPlayEnable', gf.stimDevice), 0);                 % Disable OpenEx driven sound repetition
        DA.SetTargetVal( sprintf('%s.repeatPlay',       gf.stimDevice), 0);
        
        timeNow        = DA.GetTargetVal(sprintf('%s.zTime',gf.stimDevice)) ./ gf.fStim;
        timeElapsed    = timeNow - gf.responseTime;
        timeRemaining  = gf.timeoutDuration - timeElapsed;
        
        comment = sprintf('Timeout: \nTime remaining %0.1f s', timeRemaining);
        set(h.comment,'string',comment);
        
        % Play noise
        if gf.timeoutPlayed==0
            
            % noise=rand(1,round(gf.noiseDuration/1000*gf.fStim))';
            % noise=noise .* 10^(-(30/20));
            % DA.SetTargetVal(sprintf('%s.TO_Samps', gf.stimDevice), length(noise)); % Play from central speaker
            % DA.WriteTargetVEX(sprintf('%s.speaker1', gf.stimDevice), 0, 'F32', noise'); % Play from central speaker
            
            DA.SetTargetVal(sprintf('%s.timeout',gf.stimDevice),1);
            DA.SetTargetVal(sprintf('%s.timeout',gf.stimDevice),0);
            gf.timeoutPlayed=1;
        end
        
        
        if timeRemaining <= 0,
            
            % Go straight to waiting, do not prepare another stimulus.
            % Use the same stimulus in a correction trial until the
            % correct response is made.
            % Repeat trial
            % DA.SetTargetVal( sprintf('%s.repeatPlay',       gf.stimDevice),0);
            DA.SetTargetVal( sprintf('%s.centerEnable',      gf.stimDevice), 1);
            gf.correctionTrial = 1;
            
            % Update trial number
            set(h.trialInfo,  'string',sprintf('%d', gf.TrialNumber - 1))
            gf.timeoutPlayed=0;
            gf.status  = 'WaitForStart';
        end
        
        %Update GUI
        set(h.status,'string',gf.status);
        set(h.comment,'string',comment);
        
end



%Check outputs
%checkOutputs_WE(10);                                %See toolbox for function

%Update timeline
updateTimeline_WE(20)




