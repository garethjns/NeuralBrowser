classdef SpatialMethods
    
    properties
    end
    methods (Static)
        
        function toggleSelectDataPanel(app, str)
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
            plot(app.FilteredTraceAxes, app.tvFEpoch, app.fEpoch(:,c,n))
        end
        
        function plotStim(app, n, c)
            
            % Plot on Stim tab
            % No buffering
            
            
            cla(app.StimLAAxes)
            cla(app.StimRAAxes)
            cla(app.StimLVAxes)
            cla(app.StimRVAxes)
            
            
            TT = app.Trials.TT(n,:);
            if  TT==2 || TT==4 || TT==5
                
                % Generate
                stimLA = Church2(app.stims{1,1,n});
                % Plot
                plot(app.StimLAAxes, stimLA.sound, 'Color', app.gCols(1,:))
                
                stimRA = Church2(app.stims{1,2,n});
                plot(app.StimRAAxes, stimRA.sound, 'Color', app.gCols(1,:))
                
            else
                
            end
            
            
            if TT==3 || TT==4 || TT==5
                stimLV = Church2(app.stims{2,1,n});
                plot(app.StimLVAxes, stimLV.sound, 'Color', app.gCols(2,:))
                
                stimRV = Church2(app.stims{2,2,n});
                plot(app.StimRVAxes, stimRV.sound, 'Color', app.gCols(2,:))
            else
                
            end

        end
        
        function plotLFP(app, n, c)
            switch app.GMNSwitch.Value
                case 'On'
                    norm = true;
                case 'Off'
                    norm = false;
            end
            
            plotData = app.lfpEpoch(:,c,n);
            
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
            % Calculate spikes using paramters, plot
            
            % Get channel and n
            % c = round(app.ChannelKnob.Value);
            % n = round(app.TrialSlider.Value);
            
            % Get method using selected method
            switch app.KButton.Value
                case 1
                    method = @NeuralPP.eventDetectK;
                    
                    params.medThresh = round(app.MedThreshSlider.Value);
                    params.artThresh = round(app.RejThreshSlider.Value);
                    params.eventThresh = round(app.EvThreshSlider.Value);
                    params.plotOn = false;
                    [spikes, artThresh] = method(app.fEpoch(:,c,n), ...
                        params);
                    
                    eventThresh = params.eventThresh;
                    
                    
                case 0
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
            binSize = round(app.TimeBinSlider.Value);
            
            [PSTH, ~] = NeuralAnalysis.PSTH(raster, app.fs1, binSize);
            
            % Generate tVec for this length of PSTH
            tVec = linspace(app.EpochPreTime.Value, ...
                app.EpochPostTime.Value, length(PSTH));
        end
        
        function plotPSTH(app, tVec, PSTH, n, c)
            
            % Plot PSTH
            plot(app.PSTHAxesPSTH, tVec, PSTH);
            
            % Plot stim on this tab
            SpatialMethods.plotStimOnAxes(app, app.PSTHAxesStimLeft, ...
                app.PSTHAxesStimRight, n)
            
        end
        
        function plotStimOnAxes(app, leftAxes, rightAxes, n)

            % Get times for alignmet
            preTime = app.EpochPreTime.Value;
            postTime = app.EpochPostTime.Value;
            
            cla(leftAxes)
            cla(rightAxes)
            hold(leftAxes, 'on')
            hold(rightAxes, 'on')

            TT = app.Trials.TT(n,:);
            
            leg = {};
            
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
                leg = [leg(1,:), leg(2,:)];
            end
               
            legend(leftAxes, leg)
            
        end
        
        
        function [raster, tVec] = calcRaster(app, spikes)
            [raster, tVec] = NeuralAnalysis.raster(spikes, app.fs1);
        end
        
        function plotRaster(app, raster, n, c)
            
            imagesc(app.RasterAxesRaster, double(raster));
            
            SpatialMethods.plotStimOnAxes(app, app.RasterAxesStimLeft, ...
                app.RasterAxesStimRight, n)
            
        end
        
        function app = startupFcn(app)
                       
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
                
                [fData, lfpData] = NeuralPP.filter(BB2, app.fs1);
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
                
                [fData, lfpData] = NeuralPP.filter(BB3, app.fs1);
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
        end
        
        function app = NeuralBrowseButtonPushed(app, event)
                        fn = uigetdir;
            if ischar(fn)
                app.NeuralPathEditField.Value = fn;
            end
        end
        
        function app = BehavBrowseButtonPushed(app, event)
                        [fn, pn] = uigetfile;
            if ischar(fn)
                app.BehaviouralEditField.Value = fn;
                app.behavPath = pn;
            end
        end
        
        function app = BehavLoadButtonPushed(app, event)
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
            
            for n = 1:nTrials
                disp(n)
                for v = 1:nVars
                    if isa(data{n}.(vars{v,2}), 'char')
                        trials.(vars{v,2}){n} = data{n}.(vars{v,2});
                    else
                        trials.(vars{v,2})(n,:) = data{n}.(vars{v,2});
                    end
                end
                
                % Collet stim objects to save in app.stims
                stimObjects(:,:,n) = data{n}.stimRecord;
            end
            
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
            
            fn = [pPath, block.char(), '\Epoched.mat'];
            params.EpochPreTime = app.EpochPreTime.Value;
            params.EpochPostTime = app.EpochPostTime.Value;
            if ~exist(fn, 'file')
                % If not, load processed data
                
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
                % Load existing file
                load(fn);
                app.BlockGauge.Value = 100;
            end
            
            
            app.tvFEpoch = linspace(params.EpochPreTime, ...
                params.EpochPostTime, size(fEpoch,1));
            app.tvLfpEpoch = linspace(params.EpochPreTime, ...
                params.EpochPostTime, size(lfpEpoch,1));
            
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
            
        function app = ChannelKnobValueChanged(app, event)
            value = round(app.ChannelKnob.Value);
            
            app.ChannelKnob.MajorTickLabels = {'1', '32', '63'};
            app.ChannelKnob.MajorTicks = [1, 32, 63];
            if value>1 && value<32
                app.ChannelKnob.MajorTickLabels = ...
                    {'1', num2str(value), '32', '63'};
                app.ChannelKnob.MajorTicks = [1, value, 32, 63];
            end
            
            app = SpatialMethods.ViewGroupSelectionChanged(app, event);
        end
        
        
        function app = ViewGroupSelectionChanged(app, event)
            selectedTab = app.ViewGroup.SelectedTab;
            n = round(app.TrialSlider.Value);
            c = round(app.ChannelKnob.Value);
            
            switch selectedTab.Title
                case 'Filtered trace'
                    SpatialMethods.plotFilt(app, n, c)
                case 'Spikes'
                    [spikes, artThresh, eventThresh] = ...
                        SpatialMethods.calcSpikes(app,n,c);
                    SpatialMethods.plotSpikes(app, spikes, ...
                        artThresh, eventThresh, n, c);
                    
                case 'LFP'
                    SpatialMethods.plotLFP(app, n, c);
                    
                case 'Stim'
                    SpatialMethods.plotStim(app, n, c);
                    
                case 'Raster'
                    spikes = SpatialMethods.calcSpikes(app, n, c);
                    raster = SpatialMethods.calcRaster(app, spikes);
                    SpatialMethods.plotRaster(app, raster, n, c);
                    
                case 'PSTH'
                    [spikes, artThresh, eventThresh] = ...
                        SpatialMethods.calcSpikes(app, n, c);
                    SpatialMethods.plotSpikes(app, spikes, ...
                        artThresh, eventThresh, n, c);
                    raster = SpatialMethods.calcRaster(app, spikes);
                    [PSTH, tVec] = SpatialMethods.calcPSTH(app, raster);
                    SpatialMethods.plotPSTH(app, tVec, PSTH, n,c);
                    
            end
        end
        
        function app = TrialSliderValueChanged(app, event)
            trials = app.Trials;
            
            % Get selected trial
            n = round(app.TrialSlider.Value);
            
            % Run through columns in trials and convert and non-[1, 1] 
            % to multiple columns
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
            
            app.TrialTable.Data = table2cell(trials(n,:));
            app.TrialTable.ColumnName = trials.Properties.VariableNames;
            
            SpatialMethods.ViewGroupSelectionChanged(app, event)
        end
        
    end
end