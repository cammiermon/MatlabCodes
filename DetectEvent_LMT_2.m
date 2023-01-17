function [EventDur_s,in_event,in_event_bin_nmb,in_event_bin_s,in_event_totdur,condition] = DetectEvent_LMT(fbasename,binsize,dur_session,fig,sav)

%% INPUTS 
%   fbasename     - basename of LMT behavioral events 
%                       (:,1) = Start frame of event 
%                       (:,2) = End frame of event 
%                       (:,3) = Duration of event 
%                       ex : fbasename = '3C047_SAL_NOSES' 
%                            animal_ID-condition-type event
%                       Type events 
%                           - NOSES : Nose-Nose
%                           - OGEXP : Orogenital contact by experimental mouse
%                           - OGSTM : Orogenital contact by stimulus mouse
%                           - SBSOP : Side by side contact - opposite way
%                           - APPRO : Approach by exp 
%                           - GETAW : Get away by exp 
%                           - SOESC : Social escape by exp 
%                           - SOAPP : Social approach by exp 
%                           - SBSCT : Side by side contact
%                           - BREAK : Break contact by exp 
%                           - CONTA : Contact 
%                           - GROUP : Group 2
%                       Condition
%                           - SAL : injection SALINE
%                           - OTA : injection OTA
%                           - NON : no injection (Cannula blocked one side
%                           or both sides, injection that didn't work..) 
%
%   binsize        - size of time bins to analyse data in seconds (ex : 30)
%   dur_session    - duration of recording session in minutes (ex : 10)
%   fig            - displays the figures for the session if fig ==1;
%   sav            - saves the figure and variables if sav = 1   

%% OUTPUTS 
% in_event_bin_s    - sum of the time spent in behavioral event, in second, divided in bins of 30s 
% in_event_totdur   - sum of the total time spent in behavioral event, in second
% condition         - see above for definition
% type_event        - see above for definition

%% Written by Camille Miermon 16/12/2020
%% Modified CM 09/02/21
%% Modified CM 24/04/22 add events
%% Modified CM 28/04/22 add minimum event duration, merge events if ISI too short, and add number of events 

%% Initialisation 
samplfreq = 30; %Hz
length_sess_f = dur_session*60*samplfreq; %frame

minEventDur = 0.1; %Events must be at least 0.1s long to be considered in analysis 
maxGap = 0.5; %0.5s max gap btw 2 intervals. To account for jump in detection 

animal_ID  = fbasename(1:5); 
condition = fbasename(7:9); %SAL, OTA, NON
type_event = fbasename(11:15); %OGSTM, NOSES, OGEXP, SBSOP,...
date_created = date;

%% Scrubbing data - minDur, maxGap.
dat = importdata([fbasename '.csv']);
frame = dat.data;
start_f = frame(:,1);
end_f = frame(:,2);

% Creates a vector of the size of session with 1 during event.
in_event1 = zeros(length_sess_f,1);
    for i = 1:length(start_f)   
        in_event1(start_f(i,1):end_f(i,1))=1;
    end
in_event1=in_event1(1:length_sess_f); %Put a hard threshold at dur_session
% figure, plot(in_event1), ylim([-0.5 1.5]),title('rawdata')
    
tsp = (1:1:length_sess_f)';

% Detect events
[periods,in]=Threshold([tsp,in_event1],'>',0.5);

% Scrub with maxGap 
IntNb = length(periods);
in_event2 = zeros(size(in)); 
 for j = 1:IntNb-1
     int1 = periods(j,:);
     int2 = periods(j+1,:);
     gap = int2(1,1)-int1(1,2);   
         if gap<=maxGap*samplfreq
             in_event2(int1(1,1): int2(1,2))=1; 
         else            
             in_event2(int1(1,1): int1(1,2))=1;
             in_event2(int2(1,1): int2(1,2))=1; 
         end     
 end
 
% figure, plot(in_event2), ylim([-0.5 1.5]),title('data srcub maxgap')
% [periods2,in2]=Threshold([tsp,in_event2],'>',0.5);
 
% Scrub with minDur
EventInterval = []; EventDur = []; EventNmb = 0;
L = bwlabel(in_event2);
    for ii = 1:max(L)
        idii = find(L == ii);
        dur = numel(idii);
            if dur >= minEventDur*samplfreq
                EventNmb = EventNmb+1;
                L(idii) = EventNmb;
                onset = idii(1);
                offset = idii(end);
                int = [onset offset];
                EventInterval = [EventInterval;int];
                EventDur = [EventDur;dur];
            else 
                L(idii) = 0;
            end 
    end
EventDur_s = EventDur/samplfreq;
    
% Scrubbed events:
in_event = zeros(length_sess_f,1);
    for i = 1:length(EventInterval)   
        in_event(EventInterval(i,1):EventInterval(i,2))=1;
    end  
    
% figure, plot(in_event), ylim([-0.5 1.5]),title('data srcub maxgap + minDur')
% [periods3,in3]=Threshold([tsp,in_event],'>',0.5);    
    
%% Analysis     
% Bin time of session and calculate :
    % -time spent in event 
    % -number of events 
    
binsize_f = binsize*samplfreq; 
in_event_bin = zeros(length_sess_f/binsize_f,1); %in frame
in_event_bin_nmb = zeros(length_sess_f/binsize_f,1); %in frame
startbin=1; endbin=binsize_f;
    for i = 1:length(in_event_bin)     
        in_event_bin(i,1)=sum(in_event((startbin:endbin),1));
        
        tsp =(startbin:1:endbin)';
        [periods] = Threshold([tsp,in_event(startbin:endbin)],'>',0.5);
        in_event_bin_nmb(i,1)=size(periods,1);
        
        %To avoid overestimation of event number (ie: an event falls
        % between 2 bins : 
        if i > 1
            if in_event(startbin-1) == 1
                in_event_bin_nmb(i,1) = in_event_bin_nmb(i,1)-1;
            end
        end 
        
        startbin = startbin+binsize_f;
        endbin = endbin+binsize_f;
    end
in_event_bin_s = in_event_bin/samplfreq; %in s 
in_event_totdur=sum(in_event_bin_s(:,1)); % in s  
in_event_totnmb = sum(in_event_bin_nmb);

if in_event_totnmb > EventNmb
    disp (['Overestimation of the number of events of :',num2str(in_event_totnmb-EventNmb),' events'])
elseif in_event_totnmb < EventNmb
    disp (['Underestimation of the number of events of : ',num2str(in_event_totnmb-EventNmb),' events'])
end    

if fig ==1 
    figure
    time_Axis = (binsize/60:binsize/60:dur_session)';

    subplot(2,2,1)
    %     bar(time_Axis,in_event_bin_s,'BarWidth',0.75);
    %     title (['total dur (s) =', num2str(in_event_totdur)])
    %     xlabel('Time(min)')
    %     ylabel('duration of contact (s)')
        bar(time_Axis,in_event_bin_s,'BarWidth',0.75); hold on  
        ax = gca; ax.XTick = (binsize/60*2:binsize/60*2:dur_session)'; 
        plot(time_Axis,in_event_bin_s,'color','k'), hold on
        scatter(time_Axis,in_event_bin_s,30,'k','filled')
            ax = gca; ax.XTick = (binsize/60:binsize/60:dur_session)'; 
            ylim([0 max(in_event_bin_s)+0.03])
            title (['total dur (s) =', num2str(in_event_totdur)])
            xlabel('Time(min)')
            ylabel('duration of contact (s)')

    subplot(2,2,2)
    %     plot(in_event), hold on 
    %     ylim([-0.5 1.5]),xlim([0 length_sess_f])
    %     ax = gca; ax.XTick = (binsize/60*2:binsize/60*2:dur_session)'; 
    %     title ('Events')
    %     xlabel('Time(frame)')
        plot(in_event,'lineWidth',1.5), hold on
        for j = 1:length(in_event_bin_s)
            plot([binsize_f*j binsize_f*j],[-0.25 1.25],'color','k','lineStyle','--'), hold on
        end 
            ylim([-0.5 1.5]),xlim([0 length_sess_f+1])
            ax = gca; 
            ax.XTick = (0:30*samplfreq:length_sess_f)'; 
            ax.XTickLabels = (0:0.5:10);
            xlabel('Time(min)')
            ylabel('Events')

    subplot(2,2,3)   
    %     bar(time_Axis,in_event_bin_nmb,'BarWidth',0.75);
    %     xlabel('Time(min)')
    %     ylabel('nmb of contact') 
        h = bar(time_Axis,in_event_bin_nmb,'BarWidth',0.75); hold on  
        plot(time_Axis,in_event_bin_nmb,'color','k'), hold on
        scatter(time_Axis,in_event_bin_nmb,30,'k','filled')
            h.FaceColor = 'flat'; h.CData=[0.6 0.6 0.8];
            ax = gca; ax.XTick = (binsize/60:binsize/60:dur_session)';
            ylim([0 max(in_event_bin_nmb)+1])
            title (['total nmb = ', num2str(in_event_totnmb)])
            xlabel('Time(min)')
            ylabel('number of contact')

    subplot(2,2,4) 
    % figure
        EDGES = 0:0.05:5;
        H = histogram(EventDur_s*1000,EDGES*1000);
            H.FaceColor=[1 0.4118 0.1608];
            ax = gca; 
            ax.XTick = (0:500:5000)';
            ylabel ('Number of events')
            ylim ([0 max(H.Values)+1])
            xlabel ('Duration (ms)') 
            title('Distribution of event duration')

    method = {[' minEventDur = ',num2str(minEventDur),' (s) - maxGap = ',num2str(maxGap),' (s)']}; 

    if strcmp(type_event,'NOSES') == 1 
            legend = ['Nose Nose ',num2str(animal_ID),method]; 
        elseif strcmp(type_event,'OGEXP') == 1 
             legend = ['Orogenital contact initiated by ', num2str(animal_ID),method];    
        elseif strcmp(type_event,'OGSTM') == 1 
             legend = ['Orogenital contact initiated by stim ',num2str(animal_ID),method];     
        elseif strcmp(type_event,'SBSOP') == 1
             legend = ['Side by side opposite ',num2str(animal_ID),method];  
        elseif strcmp(type_event,'APPRO') == 1
             legend = ['Approach by exp ',num2str(animal_ID),method];  
        elseif strcmp(type_event,'GETAW') == 1
             legend = ['Get away by exp ',num2str(animal_ID),method];  
        elseif strcmp(type_event,'SOESC') == 1
             legend = ['Social escape by exp ',num2str(animal_ID),method];  
        elseif strcmp(type_event,'SOAPP') == 1
             legend = ['Social approach by exp ',num2str(animal_ID),method];  
        elseif strcmp(type_event,'SBSCT') == 1
             legend = ['Side by side contact ',num2str(animal_ID),method];  
        elseif strcmp(type_event,'BREAK') == 1
             legend = ['Break contact by exp ',num2str(animal_ID),method];  
        elseif strcmp(type_event,'CONTA') == 1
             legend = ['Contact ',num2str(animal_ID),method];  
    end
    sgtitle(legend);   

end 

%% Save 
if sav == 1
    save([fbasename '_variables.mat'],'animal_ID','condition','type_event','date_created',...
        'in_event','binsize','in_event_bin','in_event_bin_s','in_event_totdur',...
        'minEventDur','maxGap','EventInterval','EventNmb','EventDur_s',...
        'in_event_bin_nmb','in_event_totnmb');
    if fig == 1 
        namefig=fbasename;
        savefig(namefig);
    end
end
end