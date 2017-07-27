classdef SpatialMethods
    % Functions from inside SpatialNeuralBrowser.mlapp. Includes callbacks.
    % Requires NeuralAnalysis library.
    % TO DO:
    % - Add batch processing functions
    % - Move out generic functions to seperate class for use with Temporal
    % interface, if created.
    
    properties
    end
    
    methods
    end
    
    methods (Static)
        
        function toggleSelectDataPanel(app, str)
            % Enable SelectDataPanel
            app.TrialTable.Enable = str;
            app.TrialSlider.Enable = str;
            app.ChannelKnob.Enable = str;
            % app.FiltSwitch.Enable = str;
            % app.PlotButton.Enable = str;
        end
        
        function updatePlots(app,n,c)
            % Clear and update plots
            % Clear all plots
            % Plot Filtered and LFP
            % (leave spikes for update button)
            % No longer used?
            
            SpatialMethods.plotFilt(app,n,c)
            SpatialMethods.plotLFP(app,n,c)
            
            [spikes, artThresh, eventThresh] = calcSpikes(app,n,c);
            SpatialMethods.plotSpikes(app, spikes, artThresh, ...
                eventThresh, n, c)
            
            SpatialMethods.plotStim(app,n,c);
            
            % And update current view
            SpatialMethods.ViewGroupSelectionChanged(app)
            
        end
        
        function plotFilt(app,n,c)
            % Plot filtered trace on appropriate axes
            plot(app.FilteredTraceAxes, app.tvFEpoch, app.fEpoch(:,c,n))
            
            % Plot stim on this tab
            SpatialMethods.plotStimOnAxes(app, ...
                app.FTStimAxes, app.FTStimAxes, n)
        end
        
        function plotStim(app, n, c)
            % Plot Stim on appropriate axes on Stim tab.
            % No epoch alignment.
            
            % Clear all axes
            cla(app.StimLAAxes)
            cla(app.StimRAAxes)
            cla(app.StimLVAxes)
            cla(app.StimRVAxes)
            
            % Get trial type
            TT = app.Trials.TT(n,:);
            
            % Plot auditory
            if  TT==2 || TT==4 || TT==5
                % Generate
                stimLA = Church2(app.stims{1,1,n});
                % Plot
                plot(app.StimLAAxes, stimLA.sound, 'Color', app.gCols(1,:))
                
                stimRA = Church2(app.stims{1,2,n});
                plot(app.StimRAAxes, stimRA.sound, 'Color', app.gCols(1,:))
                
            end
            
            % Plot visual
            if TT==3 || TT==4 || TT==5
                stimLV = Church2(app.stims{2,1,n});
                plot(app.StimLVAxes, stimLV.sound, 'Color', app.gCols(2,:))
                
                stimRV = Church2(app.stims{2,2,n});
                plot(app.StimRVAxes, stimRV.sound, 'Color', app.gCols(2,:))
            end
        end
        
        function plotLFP(app, n, c)
            % Plot LFP on LFP tab
            
            % Get global mean norm switch value
            switch app.GMNSwitch.Value
                case 'On'
                    norm = true;
                case 'Off'
                    norm = false;
            end
            
            % Get epoch
            plotData = app.lfpEpoch(:,c,n);
            
            % Normalise if on
            if norm
                plotData = plotData - int16(mean(...
                    app.lfpEpoch(:,~((1:size(app.lfpEpoch,2))==c),n), 2));
            end
            
            % Plot LFP
            plot(app.LFPAxes, app.tvLfpEpoch, plotData)
            
            % Plot stim on this tab
            SpatialMethods.plotStimOnAxes(app, ...
                app.LFPAxesStimLeft, app.LFPAxesStimRight, n)
            
        end
        
        function [spikes, artThresh, eventThresh] = calcSpikes(app,n,c)
            % Calculate spikes using paramters.
            
            % Get method using selected method
            switch app.KButton.Value
                case 1
                    method = @NeuralPP.eventDetectK;
                    
                    % Get slider values
                    params.medThresh = round(app.MedThreshSlider.Value);
                    params.artThresh = round(app.RejThreshSlider.Value);
                    params.eventThresh = round(app.EvThreshSlider.Value);
                    params.plotOn = false;
                    
                    % Calc spikes
                    [spikes, artThresh] = method(app.fEpoch(:,c,n), ...
                        params);
                    
                    eventThresh = params.eventThresh;
                    
                case 0 % 'G' - No point in using? Not implemented
                    params.thresh = [3 40];
                    params.reject = [5 50];
                    params.plotOn = false;
                    method = @NeuralPP.eventDetectG;
                    spikes = method(app.fEpoch(:,c,n), params);
                    
                    artThresh = 0;
                    eventThresh = params.thresh(3);
            end
        end
        
        function plotSpikes(app, spikes, artThresh, eventThresh, n, c)
            % Plot spikes on spikes tab
            
            hold(app.SpikesAxes, 'off')
            plot(app.SpikesAxes, app.tvFEpoch, app.fEpoch(:,c,n), ...
                'Color', [0.3, 0.3, 0.8])
            hold(app.SpikesAxes, 'on')
            plot(app.SpikesAxes, app.tvFEpoch, ...
                0.5*spikes*max(app.fEpoch(:,c,n)), ...
                'Color', [0.8, 0.1, 0.1])
            line(app.SpikesAxes, [min(app.tvFEpoch), ...
                max(app.tvFEpoch)], [eventThresh, eventThresh], ...
                'Color', 'k')
            line(app.SpikesAxes, [min(app.tvFEpoch), ...
                max(app.tvFEpoch)], [artThresh, artThresh], ...
                'Color', 'k')
            
            SpatialMethods.plotStimOnAxes(app, ...
                app.SpikesStimAxes, app.SpikesStimAxes, n)
        end
        
        function [PSTH, tVec] = calcPSTH(app, raster)
            % Calculate single-trial PSTH using specified time bin.
            
            % Get bin size
            binSize = round(app.TimeBinSlider.Value);
            
            % Calc PSTH
            [PSTH, ~] = NeuralAnalysis.PSTH(raster, app.fs1, binSize);
            
            % Generate tVec for this length of PSTH (rather than using
            % returned time vec)
            tVec = linspace(app.EpochPreTime.Value, ...
                app.EpochPostTime.Value, length(PSTH));
        end
        
        function plotPSTH(app, tVec, PSTH, n, c)
            % Plot PSTH and stim (aligned) on PSTH tab
            
            % Plot PSTH
            plot(app.PSTHAxesPSTH, tVec, PSTH);
            
            % Plot stim on this tab
            SpatialMethods.plotStimOnAxes(app, app.PSTHAxesStimLeft, ...
                app.PSTHAxesStimRight, n)
        end
        
        function plotStimOnAxes(app, leftAxes, rightAxes, n)
            % Plot stim on appropriate axes, assuming axis for each side.
            % Align stims to epoch.
            
            % Get times for alignmet
            preTime = app.EpochPreTime.Value;
            postTime = app.EpochPostTime.Value;
            
            % Clear axes
            cla(leftAxes)
            cla(rightAxes)
            hold(leftAxes, 'on')
            hold(rightAxes, 'on')
            
            % Get trial type
            TT = app.Trials.TT(n,:);
            
            % Prepare legend
            leg = {};
            
            % Prepare colours
            % If all to same axis
            % Adjust colours
            if leftAxes==rightAxes
                c1 = 1;
                c2 = 2;
                c3 = 3;
                c4 = 4;
            else
                c1 = 1;
                c2 = 1;
                c3 = 2;
                c4 = 2;
            end
            
            % Plot auditory stim
            if TT==2 || TT==4 || TT==5
                % Plot top, left stim
                % Generate
                stimLA = Church2(app.stims{1,1,n});
                [yBuff, x] = NeuralAnalysis.bufferY(stimLA.sound, ...
                    app.fs4, preTime, postTime, 0);
                
                % Plot
                plot(leftAxes, x, yBuff, 'Color', app.gCols(c1,:))
                
                % Plot bottom, right stim
                stimRA = Church2(app.stims{1,2,n});
                [yBuff, x] = NeuralAnalysis.bufferY(stimRA.sound, ...
                    app.fs4, preTime, postTime, 0);
                plot(rightAxes, x, yBuff, 'Color', app.gCols(c2,:))
                
                leg = [leg, 'Auditory'];
            end
            
            % Plot visual stim
            if TT==3 || TT==4 || TT==5
                stimLV = Church2(app.stims{2,1,n});
                [yBuff, x] = NeuralAnalysis.bufferY(stimLV.sound, ...
                    app.fs4, preTime, postTime, 0);
                plot(leftAxes, x, yBuff, 'Color', app.gCols(c3,:))
                
                stimRV = Church2(app.stims{2,2,n});
                [yBuff, x] = NeuralAnalysis.bufferY(stimRV.sound, ...
                    app.fs4, preTime, postTime, 0);
                plot(rightAxes, x, yBuff, 'Color', app.gCols(c4,:))
                leg = [leg, 'Visual'];
            end
            
            % If all to same axis
            % Adjust legend
            if leftAxes==rightAxes
                ll = 'L ' + string(leg)';
                lr = 'R ' + string(leg)';
                leg = [ll,lr];
                if size(leg,1)>1
                    leg = [leg(1,:), leg(2,:)];
                end
            end
            
            % Add legend
            legend(leftAxes, leg)
            
        end
        
        function [raster, tVec] = calcRaster(app, spikes)
            % Calc raster - single trial at the moment, and won't display
            % properly in axes
            [raster, tVec] = NeuralAnalysis.raster(spikes, app.fs1);
        end
        
        function plotRaster(app, raster, n, c)
            % Plot stim (aligned) and raster on raster tab
            % Doesn't display properly in this axes
            imagesc(app.RasterAxesRaster, double(raster));
            
            % Plot aligned stim
            SpatialMethods.plotStimOnAxes(app, app.RasterAxesStimLeft, ...
                app.RasterAxesStimRight, n)
            
        end
        
        function app = startupFcn(app)
            % Startup function
            
            % Populate available blocks in processing folder
            fn = [pwd, '\Preprocessing\*block*'];
            blocks = dir(fn);
            for b = 1:numel(blocks)
                app.BlockDropDown.Items{b} = blocks(b).name;
            end
            
        end
        
        
        % CALLBACKS
        function app = NeuralLoadButtonPushed(app, event)
            % Extract specified directory to \Extracted\ using TDTHelper
            % Process using Neural
            
            % Reset lamp and gauge
            app.NeuralLoadButton.Enable = 'off';
            app.NeuralLamp.Color = [0.9, 0.9, 0.2];
            app.NeuralGauge.Value = 0;
            
            % Try and find folder in box
            fn = app.NeuralPathEditField.Value;
            
            % Run exist checks
            if ~exist(fn, 'dir')
                app.NeuralLamp.Color = [0.9, 0, 0];
                return
            end
            
            % Try and get block name
            fps = string(fn).split('\');
            block = fps(end);
            
            if block.lower().contains('block')
                % Assume valid dir
                fn = [fn, '\'];
            else
                app.NeuralLamp.Color = [0.9, 0, 0];
                return
            end
            
            % Generate an outPath
            outPath = [pwd, '\Extracted\', block.char(), '\'];
            
            % Create TDTHelper object
            TDT = TDTHelper(fn, outPath, false);
            ok = TDT.runExtraction();
            
            app.NeuralGauge.Value = 20;
            
            if ~ok
                app.NeuralLamp.Color = [0.9, 0, 0];
                return
            end
            
            % Begin processing using NeuralPP library
            % Filter
            % Save to \Preprecessing\Block
            
            path = [pwd, '\Preprocessing\', block.char(),'\'];
            if ~exist(path, 'dir')
                mkdir(path)
            end
            
            fn = [path, 'BB2.mat'];
            if exist(fn, 'file')
                % Skip BB2
                disp('BB2 already done, skipping.')
                app.NeuralGauge.Value = 60;
            else
                BB2 = TDT.loadEvID('BB_2');
                app.NeuralGauge.Value = 30;
                
                [fData, lfpData] = ...
                    NeuralPP.filter(BB2, app.fs1, app.MemCheckBox.Value);
                app.NeuralGauge.Value = 50;
                
                save(fn, 'fData', 'lfpData')
                app.NeuralGauge.Value = 60;
            end
            clear BB2
            
            fn = [path, 'BB3.mat'];
            if exist(fn, 'file')
                % Skip BB2
                disp('BB3 already done, skipping.')
                app.NeuralGauge.Value = 100;
            else
                BB3 = TDT.loadEvID('BB_3');
                app.NeuralGauge.Value = 70;
                
                [fData, lfpData] = ...
                    NeuralPP.filter(BB3, app.fs1, app.MemCheckBox.Value);
                app.NeuralGauge.Value = 90;
                
                save(fn, 'fData', 'lfpData')
                app.NeuralGauge.Value = 100;
                
            end
            clear BB3
            
            % Update lamp
            app.NeuralLamp.Color = [0, 0.9, 0];
            
            % Add block to available block list
            if ~all(string(app.BlockDropDown.Items).contains(block))
                % If not already in list
                nItems = numel(app.BlockDropDown.Items);
                app.BlockDropDown.Items{nItems} = block.char();
            end
            
            app.NeuralLoadButton.Enable = 'on';
        end
        
        function app = NeuralBrowseButtonPushed(app, event)
            % Open file browser to find block folder
            
            fn = uigetdir;
            if ischar(fn)
                app.NeuralPathEditField.Value = fn;
            end
        end
        
        function app = BehavBrowseButtonPushed(app, event)
            % Open file brower to get behavioural file
            
            [fn, pn] = uigetfile;
            if ischar(fn)
                app.BehaviouralEditField.Value = fn;
                app.behavPath = pn;
            end
        end
        
        function app = BehavLoadButtonPushed(app, event)
            % Run on load
            
            % Reset lamp and gauge
            app.BehavLoadButton.Enable = 'off';
            app.BehaviouralLamp.Color = [0.9, 0.9, 0.2];
            app.BlockGauge.Value = 0;
            SpatialMethods.toggleSelectDataPanel(app, 'off')
            
            % Check behavioural file matches block
            block = string(app.BlockDropDown.Value);
            bFn = string(app.BehaviouralEditField.Value);
            if ~bFn.contains(block)
                app.NeuralLamp.Color = [0.9, 0, 0];
            end
            
            % Load behavioural file
            data = load([app.behavPath, bFn.char()]);
            data = data.saveData(2:end);
            app.BlockGauge.Value = 10;
            
            % Do basic processing of behav data to table
            nTrials = numel(data);
            vars = {NaN(nTrials,1), 'TT', ''; ...
                NaN(nTrials,1), 'repeatMode', ''; ...
                NaN(nTrials,1), 'side', ''; ...
                NaN(nTrials,1), 'resp', ''; ...
                NaN(nTrials,1), 'aNoise', ''; ...
                NaN(nTrials,1), 'vNoise', ''; ...
                NaN(nTrials,1), 'gap1', ''; ...
                NaN(nTrials,1), 'gap2', ''; ...
                NaN(nTrials,1), 'duration', ''; ...
                NaN(nTrials,1), 'startTrialTime', ''; ...
                NaN(nTrials,1), 'aysncOffset', 'Note typo'; ...
                cell(nTrials,1), 'eventType', ''; ...
                NaN(nTrials,2), 'speakers', ''};
            nVars = size(vars,1);
            disp(vars)
            
            trials = table(vars{:,1});
            trials.Properties.VariableNames = vars(:,2);
            trials.Properties.VariableDescriptions = vars(:,3);
            stimObjects = cell(2,2,nTrials);
            
            % For all trials
            for n = 1:nTrials
                disp(n)
                for v = 1:nVars
                    % Need to extract data from struct depending on
                    % type/shape
                    if isa(data{n}.(vars{v,2}), 'char')
                        % {}
                        trials.(vars{v,2}){n} = data{n}.(vars{v,2});
                    else
                        % ()
                        if ~isempty(data{n}.(vars{v,2}))
                            trials.(vars{v,2})(n,:) = data{n}.(vars{v,2});
                        end
                    end
                end
                
                % Collet stim objects to save in app.stims
                stimObjects(:,:,n) = data{n}.stimRecord;
                
                % Check level - gf.level unreliable
                is10 = bFn.extractBetween('level','_Rio') == '10Test';
                
                if is10 && (data{n}.TT==4 || data{n}.TT==5)
                    % There's a bug that affected saving of AV stim - vis
                    % cleared a in stimRecord. Regenerate the missing cfgs.
                    rCfg = SpatialMethods.regenCfg(data{n}, ...
                        data{n}.stimRecord);
                    
                    % Verify using vis
                    if 0
                        LV = Church2(data{n}.stimRecord{2,1});
                        RV = Church2(data{n}.stimRecord{2,2});
                        nLV = Church2(rCfg{2,1});
                        nRV = Church2(rCfg{2,2});
                        figure
                        subplot(1,2,1)
                        plot(LV.sound)
                        hold on
                        plot(RV.sound)
                        subplot(1,2,2)
                        plot(nLV.sound)
                        hold on
                        plot(nRV.sound)
                    end
                    
                    % Add this side's regrenerated A to stimRecord
                    stimObjects(1,:,n) = rCfg(1,:);
                end
            end
            
            % Save stims to app for easy access
            app.stims = stimObjects;
            clear stimObjects
            
            % Add correct column
            trials.correct = trials.side==trials.resp;
            % And save to app
            app.Trials = trials;
            
            app.BlockGauge.Value = 30;
            
            % Check if neural data already processed up to spikes
            pPath = [pwd, '\Preprocessing\'];
            if ~exist('pPath', 'dir')
                mkdir(pPath)
            end
            
            % Check if data has already been epoched
            fn = [pPath, block.char(), '\Epoched.mat'];
            params.EpochPreTime = app.EpochPreTime.Value;
            params.EpochPostTime = app.EpochPostTime.Value;
            if ~exist(fn, 'file')
                % If it hasn't, load processed data and epoch
                
                % BB2
                BB2 = load([pPath, block.char(), '\BB2.mat']);
                % Epoch
                fEpochBB2 = NeuralPP.epochData(params, ...
                    trials.startTrialTime, BB2.fData, app.fs1);
                lfpEpochBB2 = NeuralPP.epochData(params, ...
                    trials.startTrialTime, BB2.lfpData, app.fs2);
                
                app.BlockGauge.Value = 40;
                % Clean
                fEpochBB2 = NeuralPP.clean(fEpochBB2);
                app.BlockGauge.Value = 60;
                clear BB2
                
                % BB3
                BB3 = load([pPath, block.char(), '\BB3.mat']);
                
                % Epoch
                fEpochBB3 = NeuralPP.epochData(params, ...
                    trials.startTrialTime, BB3.fData, app.fs1);
                lfpEpochBB3 = NeuralPP.epochData(params, ...
                    trials.startTrialTime, BB3.lfpData, app.fs2);
                app.BlockGauge.Value = 70;
                % Clean
                fEpochBB3 = NeuralPP.clean(fEpochBB3);
                app.BlockGauge.Value = 90;
                clear BB3
                
                % Concatenate BB2 and BB3
                fEpoch = [fEpochBB2, fEpochBB3];
                clear fEpochBB2 fEpochBB3
                lfpEpoch = [lfpEpochBB2, lfpEpochBB3];
                clear lfpEpochBB2 lfpEpochBB3
                
                save(fn, 'fEpoch', 'lfpEpoch')
                app.BlockGauge.Value = 100;
                
            else
                % Load existing epoch file
                load(fn);
                app.BlockGauge.Value = 100;
            end
            
            % Save time vectors for LFP and BP data to app for easy access
            app.tvFEpoch = linspace(params.EpochPreTime, ...
                params.EpochPostTime, size(fEpoch,1));
            app.tvLfpEpoch = linspace(params.EpochPreTime, ...
                params.EpochPostTime, size(lfpEpoch,1));
            
            % Also save actual neural data to app for easy access - this
            % doesn't seem to cause noticable performance problems
            app.fEpoch = fEpoch;
            app.lfpEpoch = lfpEpoch;
            app.BehaviouralLamp.Color = [0, 0.9, 0];
            
            % Update trial viewer
            app.TrialSlider.Limits = [1, nTrials];
            SpatialMethods.TrialSliderValueChanged(app, event)
            
            % Turn on channel knob
            SpatialMethods.toggleSelectDataPanel(app, 'on')
            app.BehavLoadButton.Enable = 'on';
        end
        
        function stimRecord = regenCfg(gf, stimRecord)
            % Regenerate missing auditory cfgs from incomplete stimRecord
            % and gf snapshot
            % Get:
            % - common settings from gf
            % - Side specific settings from visual stims in stimRecord
            %
            % Retruns complete stimRecord with regenerated V stimuli as
            % well. These can/should be used for verification in calling
            % function!
            %
            % Undoes this bug in level 10:
            %
            %-------------------
            % .... [Generate] ....
            % if isstruct(aStim)
            %    soundLR(1:length(aStim.sound),i) = ...
            %        aStim.sound';
            %    aStim = rmfield(aStim, 'sound');
            %     gf.stimRecord{1,i}=aStim;
            %     gf.stimRecord{2,i}=NaN;     <---- Clear previous
            % end
            % if isstruct(vStim)
            %     lightLR(1:length(vStim.sound),i) = ...
            %         (vStim.sound*gf.vMulti)';
            %     vStim = rmfield(vStim, 'sound');
            %     gf.stimRecord{2,i}=vStim;
            %     gf.stimRecord{1,i}=NaN;    <---- Clear A just set
            % end
            % --------------------
            
            % Therefore, (may) need to recreate A stim cfg while
            % loading
            % Relevant info should be in gf, snapshot saved in trial
            % information
            
            % Missing gf.rideNoise?
            % Get from other stim objects
            gf.rideNoise = stimRecord{2,1}.rideNoise;
            % nEvents also not saved in gf, get from corresponding
            % vis for each side:
            eV = [stimRecord{2,1}.nEvents, ...
                stimRecord{2,2}.nEvents];
            
            for s = 1:2
                
                % Common parameters
                % p=gf.speakers(i);
                % (here left=8 and right=2, so ev processed in correct order)
                cfg.nEvents = eV(s);
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
                
                
                % Using code from level10Test_WE.m
                switch gf.TT
                    case 4 % Sync
                        cfg.type = 'Aud';
                        cfg.rise = gf.stimRise;
                        cfg.engBuff = gf.endBuff;
                        cfg.startBuff = gf.startBuff;
                        cfg.noiseMag = gf.aNoise;
                        cfg.rideNoise = gf.rideNoise;
                        % Seed
                        cfg.seed = gf.seeds(3,1);
                        
                        aStim = Church2(cfg);
                        
                        % Keep cfg, replace V fields
                        cfg.type = 'Vis';
                        cfg.eventType = 'flat';
                        cfg.noiseMag = gf.vNoise;
                        cfg.eventMag = gf.vRange;
                        % Seed
                        cfg.seed = gf.seeds(3,2); % (Same as A)
                        
                        %
                        vStim = Church2(cfg);
                    case 5 % Async
                        
                        cfg.type = 'Aud';
                        cfg.engBuff = gf.endBuff;
                        cfg.startBuff = gf.startBuff;
                        cfg.noiseMag = gf.aNoise;
                        cfg.rideNoise = gf.rideNoise;
                        % Seed
                        cfg.seed = gf.seeds(4,1);
                        
                        aStim = Church2(cfg);
                        
                        clear cfg
                        
                        % Re-set common
                        cfg.Fs = gf.fStim; % Sampling rate
                        cfg.eventLength = gf.stimEventDuration;
                        % cfg.mag = 0-gf.atten;
                        cfg.eventType = gf.eventType;
                        cfg.gap1 = gf.gap1; % Gap durations
                        cfg.gap2 = gf.gap2;
                        cfg.duration=gf.duration;
                        cfg.cutOff=gf.cutOff;
                        cfg.eventFreq = ...
                            gf.eventFreq; % If 'sine', frequency
                        cfg.noiseType='multipleBlocks';
                        cfg.rise=gf.stimRise;
                        cfg.rideNoise=gf.rideNoise;
                        
                        % Disabled - don't want to reset in regen code,
                        % rather just use existing value
                        % Select variable offset
                        % gf.aysncOffset = ...
                        %     gf.asyncOffsets(randi(numel(...
                        %     gf.asyncOffsets)));
                        
                        % Set specific
                        cfg.type = 'Vis';
                        cfg.eventType = 'flat';
                        cfg.eventMag = gf.vRange;      % From above
                        cfg.endBuff = gf.endBuff + gf.aysncOffset;
                        cfg.startBuff = gf.startBuff + gf.aysncOffset;
                        cfg.noiseMag = gf.vNoise;
                        % Seed
                        cfg.seed = gf.seeds(4,2);
                        
                        vStim=Church2(cfg);
                end
                
                stimRecord{1,s} = rmfield(aStim, 'sound');
                stimRecord{2,s} = rmfield(vStim, 'sound');
            end
            
        end
        
        function app = ChannelKnobValueChanged(app, event)
            % Run on channel selection change
            
            % Get value
            value = round(app.ChannelKnob.Value);
            
            % Set major tick labels to this combination to display "1" and
            % "32" correctly
            app.ChannelKnob.MajorTickLabels = {'1', '32', '63'};
            app.ChannelKnob.MajorTicks = [1, 32, 63];
            
            % Then add the current selection to the tick labels so choice
            % is obvious
            if value>1 && value<32
                app.ChannelKnob.MajorTickLabels = ...
                    {'1', num2str(value), '32', '63'};
                app.ChannelKnob.MajorTicks = [1, value, 32, 63];
            end
            
            % Recalc and replot current tab
            app = SpatialMethods.ViewGroupSelectionChanged(app, event);
        end
        
        function app = ViewGroupSelectionChanged(app, event)
            % This function relacs and replots for currentely selected view
            % tab
            
            % Get tab, trial, channel info
            selectedTab = app.ViewGroup.SelectedTab;
            n = round(app.TrialSlider.Value);
            c = round(app.ChannelKnob.Value);
            
            % Switch of slected tab to limit processing to current only
            switch selectedTab.Title
                case 'Filtered trace'
                    % Plot the filtered trace
                    SpatialMethods.plotFilt(app, n, c)
                    
                case 'Spikes'
                    % Calc and plot spikes
                    [spikes, artThresh, eventThresh] = ...
                        SpatialMethods.calcSpikes(app,n,c);
                    SpatialMethods.plotSpikes(app, spikes, ...
                        artThresh, eventThresh, n, c);
                    
                case 'LFP'
                    % Plot LFP
                    SpatialMethods.plotLFP(app, n, c);
                    
                case 'Stim'
                    % Plot unaliged stim on seperate axes
                    SpatialMethods.plotStim(app, n, c);
                    
                case 'Raster'
                    % Calc spikes, calc and plot the pointless raster
                    spikes = SpatialMethods.calcSpikes(app, n, c);
                    raster = SpatialMethods.calcRaster(app, spikes);
                    SpatialMethods.plotRaster(app, raster, n, c);
                    
                case 'PSTH'
                    % Calc spikes, calc raster, calc and plot PSTH
                    [spikes, artThresh, eventThresh] = ...
                        SpatialMethods.calcSpikes(app, n, c);
                    % And plot spikes??
                    SpatialMethods.plotSpikes(app, spikes, ...
                        artThresh, eventThresh, n, c);
                    
                    raster = SpatialMethods.calcRaster(app, spikes);
                    
                    [PSTH, tVec] = SpatialMethods.calcPSTH(app, raster);
                    SpatialMethods.plotPSTH(app, tVec, PSTH, n,c);
            end
        end
        
        function app = TrialSliderValueChanged(app, event)
            % Run on trial selction change
            
            trials = app.Trials;
            
            % Get selected trial
            n = round(app.TrialSlider.Value);
            
            % Run through columns in trials and convert and non-[1, 1]
            % to multiple columns
            % Drop anything that can't be displayed in UI element
            m = width(trials);
            dropList = {};
            newTable = table();
            for c = 1:m
                m2 = size(trials{:,c},2);
                if m2>1
                    dropList = [dropList, ...
                        trials.Properties.VariableNames{c}];
                    for c2 = 1:m2
                        newName = [trials.Properties.VariableNames{c}, ...
                            '_', num2str(c2)];
                        
                        newTable.(newName) = trials{:,c}(:,c2);
                    end
                end
                
            end
            
            trials(:,dropList) = [];
            trials = [trials, newTable];
            
            % Update table
            app.TrialTable.Data = table2cell(trials(n,:));
            app.TrialTable.ColumnName = trials.Properties.VariableNames;
            
            % And clac/plot current view tab
            SpatialMethods.ViewGroupSelectionChanged(app, event)
        end
        

    end
end