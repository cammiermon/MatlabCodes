function [in_Contact_OTA,in_Contact_SAL,in_Contact_SAL_nmb,in_Contact_OTA_nmb]= AnalyseEvent_LMT_2(binsize,dur_session,fig,sav)

%% INPUTS 
%   binsize        - size of time bins to analyse data in seconds (ex : 30)
%   dur_session    - duration of recording session in minutes (ex : 10)
%   fig            - displays the figures if fig ==1;
%   sav            - saves the figure and variables if sav = 1  
%
%  This code:
%   - is written to analyse LMT data of mice who received OTA, SAL, or NON injection 
%   - should be run in a folder with data of ONLY one behavioral event : 
%     * NOSES (Nose-Nose) or
%     * OGEXP (Orogenital from exp) or
%     * OGSTM (Orogenital from stim) or
%     * SBSOP (Side by side opposite) or 
%     * APPRO (Approach by exp) or 
%     * GETAW (Get away by exp) or
%     * SOESC (Social escape by exp) or 
%     * SOAPP (Social approach by exp) or
%     * SBSCT (Side by side contact) or 
%     * BREAK (Break contact by exp) or 
%     * CONTA (Contact) or
%     * GROUP (Group 2)
%   - saves individual data/figures if the sav variable of DetectEvent_LMT
%     is equal to 1;
%
%% DEPENDENCY 
% DetectEvent_LMT
% or DetectEvent_LMT_2 (with scrubbed data - minDur, maxGap, and number of events
%
%% Written by Camille Miermon 16/12/2020
%% Modified by CM with GC inputs 16/02/2021
%% Modified by CM add more events 24/04/2022

%% Initialization 

if ~exist('binsize','var') ; binsize = 30; else ; if isempty(binsize) ; binsize =30; end ; end
if ~exist('dur_session','var') ; dur_session = 10; else ; if isempty(dur_session) ; dur_session = 10; end ; end
if ~exist('fig','var') ; fig = 1; else ; if isempty(fig) ; fig =1; end ; end
if ~exist('sav','var') ; sav = 0; else ; if isempty(sav) ; sav =0; end ;end

date_created = date;

t_BinSize = binsize/60; %min
timeVector = transpose(t_BinSize/2:t_BinSize:dur_session);
t_BinSizeIndices = transpose(1:1:dur_session/t_BinSize);

%% Get all sessions 
myDir = uigetdir; %gets directory
cd(myDir);
SessionAnalysed = dir(fullfile(myDir,'*.csv')); %gets all csv files in struct

in_Contact = []; in_Contact_totdur=[]; in_Contact_nmb = []; in_Contact_totnmb=[];
Dur_OTA = []; Dur_SAL = []; Dur_NON = [];
GroupID = []; GroupID2 = []; timestamps=[]; i=1;j=1;l=1;
in_Contact_SAL =[];in_Contact_SAL_nmb=[];
in_Contact_OTA =[];in_Contact_OTA_nmb=[];
in_Contact_NON =[];in_Contact_NON_nmb=[];

for k = 1:length(SessionAnalysed)
  baseFileName = SessionAnalysed(k).name;
  fprintf(1, 'Now reading %s\n', baseFileName);
  fbasename = baseFileName(1:15);

  [EventDur_s,in_event,in_event_bin_nmb,in_event_bin_s,in_event_totdur,condition]...
      = DetectEvent_LMT_2(fbasename,binsize,dur_session,0,sav);
 
  in_Contact = [in_Contact;in_event_bin_s];
  in_Contact_totdur = [in_Contact_totdur;in_event_totdur];
  in_Contact_nmb = [in_Contact_nmb;in_event_bin_nmb];
  in_Contact_totnmb = [in_Contact_totnmb;sum(in_event_bin_nmb)];
    if strcmp(condition,'OTA') == 1 
          Group =  1;
          in_Contact_OTA(:,i)= in_event_bin_s;  
          in_Contact_OTA_nmb(:,i)= in_event_bin_nmb;
          Dur_OTA=[Dur_OTA;EventDur_s];
          i=i+1;
    elseif strcmp(condition,'SAL') == 1
          Group =  2;
          in_Contact_SAL(:,j)= in_event_bin_s;  
          in_Contact_SAL_nmb(:,j)= in_event_bin_nmb;
          Dur_SAL=[Dur_SAL;EventDur_s];
          j=j+1; 
    elseif strcmp(condition,'NON') == 1
          Group =  3;
          in_Contact_NON(:,l)= in_event_bin_s; 
          in_Contact_NON_nmb(:,l)= in_event_bin_nmb;
          Dur_NON=[Dur_NON;EventDur_s];
          l=l+1;
    end 
  GroupID = [GroupID ; repmat(Group,numel(in_event_bin_s),1)];
  GroupID2 = [GroupID2; Group];
  timestamps = [timestamps;t_BinSizeIndices];   
end

type_event = fbasename(11:15);
    if strcmp(type_event,'NOSES') == 1 
        Type_event = ('Nose Nose');
    elseif strcmp(type_event,'OGEXP') == 1
        Type_event = ('Orogenital contact initiated by exp');   
    elseif strcmp(type_event,'OGSTM') == 1 
        Type_event = ('Orogenital contact initiated by stim');  
    elseif strcmp(type_event,'SBSOP') == 1
        Type_event = ('Side by side opposite');
    elseif strcmp(type_event,'APPRO') == 1
        Type_event = ('Approach by exp');
    elseif strcmp(type_event,'GETAW') == 1
        Type_event = ('Get away by exp');
    elseif strcmp(type_event,'SOESC') == 1
        Type_event = ('Social escape by exp');
    elseif strcmp(type_event,'SOAPP') == 1
        Type_event = ('Social approach by exp');
    elseif strcmp(type_event,'SBSCT') == 1
        Type_event = ('Side by side contact'); 
    elseif strcmp(type_event,'BREAK') == 1
        Type_event = ('Break contact by exp'); 
    elseif strcmp(type_event,'CONTA') == 1
        Type_event = ('Contact');
    elseif strcmp(type_event,'GROUP') == 1
        Type_event = ('Group 2'); 
    end
    
%% Analysis individual data 
in_Contact = [in_Contact GroupID timestamps];
in_Contact_totdur = [in_Contact_totdur GroupID2];
in_Contact_nmb = [in_Contact_nmb GroupID timestamps];
in_Contact_totnmb = [in_Contact_totnmb GroupID2];

in_Contact_OTA_cum = cumsum(in_Contact_OTA); % Cumulative data
in_Contact_SAL_cum = cumsum(in_Contact_SAL);
    if ~isempty(in_Contact_NON)
        in_Contact_NON_cum = cumsum(in_Contact_NON);
    end    

in_Contact_OTA_nmb_cum = cumsum(in_Contact_OTA_nmb);
in_Contact_SAL_nmb_cum = cumsum(in_Contact_SAL_nmb);
     if ~isempty(in_Contact_NON_nmb)
        in_Contact_NON_nmb_cum = cumsum(in_Contact_NON_nmb);
    end  

%% Analysis grouped data
Group_N = [i-1 j-1 l-1]; % n number of animals for OTA, SAL, NON
Group_N = Group_N(Group_N>0);

% Duration data 
MeanAcrossGroups = transpose(accumarray(in_Contact(:,2:3) , in_Contact(:,1) , [] , @nanmean)); %Bin data
StdAcrossGroups = transpose(accumarray(in_Contact(:,2:3) , in_Contact(:,1) , [] , @nanstd));
SEMAcrossGroups = StdAcrossGroups./sqrt(Group_N);

MeanAcrossGroups_cum = cumsum(MeanAcrossGroups); %Cumulative data
in_Contact_OTA_cum_SEM = std(in_Contact_OTA_cum,[],2)./sqrt(Group_N(:,1));  
in_Contact_SAL_cum_SEM = std(in_Contact_SAL_cum,[],2)./sqrt(Group_N(:,2));
    if ~isempty(in_Contact_NON)
        in_Contact_NON_cum_SEM = std(in_Contact_NON_cum,[],2)./sqrt(Group_N(:,3));
    end

MeanAcrossGroups_totdur = transpose(accumarray(in_Contact_totdur(:,2) , in_Contact_totdur(:,1) , [] , @nanmean));
StdAcrossGroups_totdur = transpose(accumarray(in_Contact_totdur(:,2) , in_Contact_totdur(:,1) , [] , @nanstd));
SEMAcrossGroups_totdur = StdAcrossGroups_totdur./sqrt(Group_N);

% Number of events data 
MeanAcrossGroups_2 = transpose(accumarray(in_Contact_nmb(:,2:3) , in_Contact_nmb(:,1) , [] , @nanmean)); %Bin data
StdAcrossGroups_2 = transpose(accumarray(in_Contact_nmb(:,2:3) , in_Contact_nmb(:,1) , [] , @nanstd));
SEMAcrossGroups_2 = StdAcrossGroups_2./sqrt(Group_N);

MeanAcrossGroups_cum_2 = cumsum(MeanAcrossGroups_2); %Cumulative data
in_Contact_OTA_cum_SEM_2 = std(in_Contact_OTA_nmb_cum,[],2)./sqrt(Group_N(:,1));  
in_Contact_SAL_cum_SEM_2 = std(in_Contact_SAL_nmb_cum,[],2)./sqrt(Group_N(:,2));
    if ~isempty (in_Contact_NON)
        in_Contact_NON_cum_SEM_2 = std(in_Contact_NON_nmb_cum,[],2)./sqrt(Group_N(:,3));
    end

MeanAcrossGroups_totnmb = transpose(accumarray(in_Contact_totnmb(:,2) , in_Contact_totnmb(:,1) , [] , @nanmean));
StdAcrossGroups_totnmb = transpose(accumarray(in_Contact_totnmb(:,2) , in_Contact_totnmb(:,1) , [] , @nanstd));
SEMAcrossGroups_totnmb = StdAcrossGroups_totnmb./sqrt(Group_N);


%% figures 

if fig ==1
    
close all
    %% figure individual data     
%Time spent in behavioral event per time bin
figure(1); 
    plot (timeVector,in_Contact_OTA,'Color','r'), hold on
    plot (timeVector,in_Contact_SAL,'Color','k'), hold on
    if ~isempty(in_Contact_NON)
        plot (timeVector,in_Contact_NON,'Color',[0.5,0.5,0.5]), hold off
    end
    xlabel('Time(min)')
    ylabel('duration of contact (s)')
    title([num2str(Type_event),'- individual data dur'])
    
namefig1=[num2str(Type_event),'- individual data'];  
    
%Cumulative time spent in behavioral event
figure(2); 
    plot (timeVector,in_Contact_OTA_cum,'Color','r'), hold on
    plot (timeVector,in_Contact_SAL_cum,'Color','k'), hold on
    if ~isempty(in_Contact_NON)
        plot (timeVector,in_Contact_NON_cum,'Color',[0.5,0.5,0.5]), hold off
    end
    xlabel('Time(min)')
    ylabel('cumulative time spent in contact (s)')
    title([num2str(Type_event),'- individual data']) 
    
namefig2=[num2str(Type_event),'- individual data cumu dur'];

% Figure grouped data 

%Mean data per group 
figure(3);
    errorbar(timeVector,MeanAcrossGroups(:,1),SEMAcrossGroups(:,1),'Color','r'), hold on 
    errorbar(timeVector,MeanAcrossGroups(:,2),SEMAcrossGroups(:,2),'Color','k'), hold on 
    if length(Group_N)>2
        errorbar(timeVector,MeanAcrossGroups(:,3),SEMAcrossGroups(:,3),'Color',[0.5,0.5,0.5])
    end
    xlabel('Time(min)')
    ylabel('duration of contact (s)')
    title([num2str(Type_event),'- Mean data +/- SEM'])    

namefig3=[num2str(Type_event),'- Mean data dur'];

%Cumulative Mean data per group     
figure(4);
    errorbar(timeVector,MeanAcrossGroups_cum(:,1),in_Contact_OTA_cum_SEM,'Color','r'), hold on 
    errorbar(timeVector,MeanAcrossGroups_cum(:,2),in_Contact_SAL_cum_SEM,'Color','k'), hold on 
    if length(Group_N)>2
        errorbar(timeVector,MeanAcrossGroups_cum(:,3),in_Contact_NON_cum_SEM,'Color',[0.5,0.5,0.5])
    end
    xlabel('Time(min)')
    ylabel('cumulative time spent in contact (s)')
    title([num2str(Type_event),'- Mean data'])    

namefig4=[num2str(Type_event),'- Mean data cumu dur']; 

%Mean data total time per group 
figure(5); 
    y_axis=(1:1:numel(MeanAcrossGroups_totdur));
    a = bar(y_axis,MeanAcrossGroups_totdur,'BarWidth',0.75); hold on 
    errorbar(y_axis,MeanAcrossGroups_totdur,SEMAcrossGroups_totdur, 'k', 'linestyle', 'none'); %SEM
    a.FaceColor = 'flat';
    a.CData(1,:)=[1,0,0];
    a.CData(2,:)=[0,0,0];
    if length(Group_N)>2
        a.CData(3,:)=[0.5,0.5,0.5];
    end
    ax = gca;  
    if length(Group_N)>2
        ax.XTick = [1 2]; 
        ax.XTickLabels = {'OTA','SAL'};
    else
        ax.XTick = [1 2 3]; 
        ax.XTickLabels = {'OTA','SAL','NON'};
       
    end
    ylim([0,max(MeanAcrossGroups_totdur)+4]);
    xlim([0.4 numel(MeanAcrossGroups_totdur)+0.6])
    ylabel('Total dur of contact (s)');
    title([num2str(Type_event),'- Mean data']) 
    
namefig5=[num2str(Type_event),'- Mean data barplot - dur']; 

% Summary main data 
figure(6);
    subplot(1,3,1)
        plot (timeVector,in_Contact_OTA_cum,'Color','r'), hold on
        plot (timeVector,in_Contact_SAL_cum,'Color','k'), hold on
        if ~isempty(in_Contact_NON)
            plot (timeVector,in_Contact_NON_cum,'Color',[0.5,0.5,0.5]), hold on
        end
        xlabel('Time(min)')
        ylabel('cumulative time spent in contact (s)')
        title('Individual cumulative data') 
    subplot(1,3,2)
        errorbar(timeVector,MeanAcrossGroups_cum(:,1),in_Contact_OTA_cum_SEM,'Color','r'), hold on 
        errorbar(timeVector,MeanAcrossGroups_cum(:,2),in_Contact_SAL_cum_SEM,'Color','k'), hold on 
        if ~isempty(in_Contact_NON)
        	errorbar(timeVector,MeanAcrossGroups_cum(:,3),in_Contact_NON_cum_SEM,'Color',[0.5,0.5,0.5])
        end
        xlabel('Time(min)')
        ylabel('cumulative time spent in contact (s)')
        title('Mean cumululative data +/- SEM')
    subplot(1,3,3)
        y_axis=(1:1:numel(MeanAcrossGroups_totdur));
        a = bar(y_axis,MeanAcrossGroups_totdur,'BarWidth',0.75); hold on 
        errorbar(y_axis,MeanAcrossGroups_totdur,SEMAcrossGroups_totdur, 'k', 'linestyle', 'none'); %SEM
        a.FaceColor = 'flat';
        a.CData(1,:)=[1,0,0];
        a.CData(2,:)=[0,0,0];
        if ~isempty(in_Contact_NON)
            a.CData(3,:)=[0.5,0.5,0.5];
        end 
        ax = gca; 
        if length(Group_N)>2
            ax.XTick = [1 2 3]; 
            ax.XTickLabels = {'OTA','SAL','NON'};
        else 
            ax.XTick = [1 2]; 
            ax.XTickLabels = {'OTA','SAL'};
        end
        ylim([0,max(MeanAcrossGroups_totdur)+4]);
        xlim([0.4 numel(MeanAcrossGroups_totdur)+0.6])
        ylabel('Total dur of contact (s)');
        title('Total time spent +/- SEM')
        
namefig6=[num2str(Type_event),'- main plots - dur'];

if ~isempty(in_Contact_NON)
    legend = [num2str(Type_event),'; OTA = ',num2str(Group_N(:,1)),', SAL =',num2str(Group_N(:,2)),', NON =',num2str(Group_N(:,3))]; 
else
    legend = [num2str(Type_event),'; OTA = ',num2str(Group_N(:,1)),', SAL =',num2str(Group_N(:,2))];
end
sgtitle(legend);
    

% Summary main data + nmb + distribution of event duration
figure(7); 
     subplot(3,3,1)
        plot (timeVector,in_Contact_OTA_cum,'Color','r'), hold on
        plot (timeVector,in_Contact_SAL_cum,'Color','k'), hold on
        if ~isempty(in_Contact_NON)
            plot (timeVector,in_Contact_NON_cum,'Color',[0.5,0.5,0.5]), hold on
        end
        xlabel('Time(min)')
        ylabel('Cumulative time (s)')
        title('Individual data') 
    subplot(3,3,4)
        errorbar(timeVector,MeanAcrossGroups_cum(:,1),in_Contact_OTA_cum_SEM,'Color','r'), hold on 
        errorbar(timeVector,MeanAcrossGroups_cum(:,2),in_Contact_SAL_cum_SEM,'Color','k'), hold on 
        if ~isempty(in_Contact_NON)
            errorbar(timeVector,MeanAcrossGroups_cum(:,3),in_Contact_NON_cum_SEM,'Color',[0.5,0.5,0.5])
        end
        xlabel('Time(min)')
        ylabel('Cumulative time (s)')
        title('Group data')
    subplot(3,3,7)
        y_axis=(1:1:numel(MeanAcrossGroups_totdur));
        a = bar(y_axis,MeanAcrossGroups_totdur,'BarWidth',0.75); hold on 
        errorbar(y_axis,MeanAcrossGroups_totdur,SEMAcrossGroups_totdur, 'k', 'linestyle', 'none'); %SEM
        a.FaceColor = 'flat';
        a.CData(1,:)=[1,0,0];
        a.CData(2,:)=[0,0,0];
        if ~isempty(in_Contact_NON)
            a.CData(3,:)=[0.5,0.5,0.5];
        end
        ax = gca; 
        if length(Group_N)>2
            ax.XTick = [1 2 3]; 
            ax.XTickLabels = {'OTA','SAL','NON'};
        else 
            ax.XTick = [1 2]; 
            ax.XTickLabels = {'OTA','SAL'};
        end
        ylim([0,max(MeanAcrossGroups_totdur)+(max(MeanAcrossGroups_totdur)-1)]);
        xlim([0.4 numel(MeanAcrossGroups_totdur)+0.6])
        ylabel('Total dur of contact (s)');
        title('Group data')
    subplot(3,3,2)
        plot (timeVector,in_Contact_OTA_nmb_cum,'Color','r'), hold on
        plot (timeVector,in_Contact_SAL_nmb_cum,'Color','k'), hold on
        if ~isempty(in_Contact_NON)
            plot (timeVector,in_Contact_NON_nmb_cum,'Color',[0.5,0.5,0.5]), hold on
        end
        xlabel('Time(min)')
        ylabel('Cumulative number of event')
        title('Individual data') 
    subplot(3,3,5)
        errorbar(timeVector,MeanAcrossGroups_cum_2(:,1),in_Contact_OTA_cum_SEM_2,'Color','r'), hold on 
        errorbar(timeVector,MeanAcrossGroups_cum_2(:,2),in_Contact_SAL_cum_SEM_2,'Color','k'), hold on 
        if ~isempty(in_Contact_NON)
            errorbar(timeVector,MeanAcrossGroups_cum_2(:,3),in_Contact_NON_cum_SEM_2,'Color',[0.5,0.5,0.5])
        end
        xlabel('Time(min)')
        ylabel('Cumulative number of event')
        title('Group data')
    subplot(3,3,8)
        y_axis=(1:1:numel(MeanAcrossGroups_totnmb));
        a = bar(y_axis,MeanAcrossGroups_totnmb,'BarWidth',0.75); hold on 
        errorbar(y_axis,MeanAcrossGroups_totnmb,SEMAcrossGroups_totnmb, 'k', 'linestyle', 'none'); %SEM
        a.FaceColor = 'flat';
        a.CData(1,:)=[1,0,0];
        a.CData(2,:)=[0,0,0];
        if ~isempty(in_Contact_NON)
            a.CData(3,:)=[0.5,0.5,0.5];
        end
        ax = gca;
        if length(Group_N)>2
            ax.XTick = [1 2 3]; 
            ax.XTickLabels = {'OTA','SAL','NON'};
        else 
           ax.XTick = [1 2]; 
           ax.XTickLabels = {'OTA','SAL'};
        end
        ylim([0,max(MeanAcrossGroups_totnmb)+(max(SEMAcrossGroups_2,[],'all')+6)]);
        xlim([0.4 numel(MeanAcrossGroups_totnmb)+0.6])
        ylabel('Total nmb of events');
        title('Group data')   

    EDGES = 0:0.05:5;
    subplot(3,3,3)
    H = histogram(Dur_OTA*1000,EDGES*1000,'Normalization','probability');
        H.FaceColor=[1,0,0];
        ax = gca; 
        ax.XTick = (0:500:5000)';
        ylabel ('Normalised number of events')
        ylim([0 0.2])
        xlabel ('Duration (ms)') 
        title('Distribution of event duration - OTA')
    subplot(3,3,6)
    H = histogram(Dur_SAL*1000,EDGES*1000,'Normalization','probability');
        H.FaceColor=[0,0,0];
        ax = gca; 
        ax.XTick = (0:500:5000)';
        ylabel ('Normalised number of events')
        ylim([0 0.2])
        xlabel ('Duration (ms)') 
        title('Distribution of event duration - SAL')
    subplot(3,3,9)
    H = histogram(Dur_NON*1000,EDGES*1000,'Normalization','probability');
        H.FaceColor=[0.5,0.5,0.5];
        ax = gca; 
        ax.XTick = (0:500:5000)';
        ylabel ('Normalised number of events')
        ylim([0 0.2])
        xlabel ('Duration (ms)') 
        title('Distribution of event duration - NON')
        
namefig7=[num2str(Type_event),'- main plots - dur - nmb'];

if length(Group_N)>2
    legend = [num2str(Type_event),' OTA = ',num2str(Group_N(:,1)),', SAL =',num2str(Group_N(:,2)),', NON =',num2str(Group_N(:,3))];
else 
    legend = [num2str(Type_event),' OTA = ',num2str(Group_N(:,1)),', SAL =',num2str(Group_N(:,2))];
end 
sgtitle(legend);    

% Distribution of event duration 
    EDGES = 0:0.05:5;
    figure(20),A1 = histogram(Dur_OTA*1000,EDGES*1000,'Normalization','probability');
     figure(21),A2 = histogram(Dur_SAL*1000,EDGES*1000,'Normalization','probability');
      figure(22),A3 = histogram(Dur_NON*1000,EDGES*1000,'Normalization','probability');
      
    figure(8)
        plot(EDGES(2:end)*1000,A1.Values,'Color','r'), hold on 
        plot(EDGES(2:end)*1000,A2.Values,'Color','k'), hold on 
        plot(EDGES(2:end)*1000,A3.Values,'Color',[0.5,0.5,0.5]) 
            ylabel ('Normalised number of events')
            xlim([0 3000]), ylim([0 0.2])
            xlabel ('Duration (ms)') 
            if length(Group_N)>2
                legend = [num2str(Type_event),' OTA = ',num2str(Group_N(:,1)),', SAL =',num2str(Group_N(:,2)),', NON =',num2str(Group_N(:,3))];
            else 
                legend = [num2str(Type_event),' OTA = ',num2str(Group_N(:,1)),', SAL =',num2str(Group_N(:,2))];
            end
            sgtitle(legend); 
            namefig8=[num2str(Type_event),'- distribution event dur '];

                close(figure(20)), close(figure(21)), close(figure(22))
    
% Individual data - number of event per time bin 
figure(9); 
    errorbar(timeVector,MeanAcrossGroups_2(:,1),SEMAcrossGroups_2(:,1),'Color','r'), hold on 
    errorbar(timeVector,MeanAcrossGroups_2(:,2),SEMAcrossGroups_2(:,2),'Color','k'), hold on 
    if ~isempty(in_Contact_NON)
        errorbar(timeVector,MeanAcrossGroups_2(:,3),SEMAcrossGroups_2(:,3),'Color',[0.5,0.5,0.5])
    end
    xlabel('Time(min)')
    ylabel('number of events')
    title([num2str(Type_event),' - number of events'])  
    
namefig9=[num2str(Type_event),'- mean nmb'];    

end

%% Sav
if sav ==1  
    if ~isempty(in_Contact_NON)
    save([num2str(Type_event) '_variables.mat'],'type_event','timeVector','date_created','binsize','Group_N',...
        'in_Contact','in_Contact_NON','in_Contact_OTA','in_Contact_SAL','in_Contact_OTA_cum','in_Contact_SAL_cum',...
        'in_Contact_NON_cum','in_Contact_totdur','MeanAcrossGroups','MeanAcrossGroups_cum','MeanAcrossGroups_totdur',...
        'SessionAnalysed','in_Contact_nmb','in_Contact_NON_nmb','in_Contact_SAL_nmb','in_Contact_OTA_nmb',...
        'in_Contact_totnmb','in_Contact_OTA_nmb_cum','in_Contact_SAL_nmb_cum','in_Contact_NON_nmb_cum',...
        'Dur_NON','Dur_OTA','Dur_SAL');
    else 
        save([num2str(Type_event) '_variables.mat'],'type_event','timeVector','date_created','binsize','Group_N',...
        'in_Contact','in_Contact_OTA','in_Contact_SAL','in_Contact_OTA_cum','in_Contact_SAL_cum',...
        'in_Contact_totdur','MeanAcrossGroups','MeanAcrossGroups_cum','MeanAcrossGroups_totdur',...
        'SessionAnalysed','in_Contact_nmb','in_Contact_SAL_nmb','in_Contact_OTA_nmb',...
        'in_Contact_totnmb','in_Contact_OTA_nmb_cum','in_Contact_SAL_nmb_cum',...
        'Dur_OTA','Dur_SAL');
    end 
    if fig ==1
        savefig(figure(1),namefig1);
        savefig(figure(2),namefig2);
        savefig(figure(3),namefig3);
        savefig(figure(4),namefig4);
        savefig(figure(5),namefig5);
        savefig(figure(6),namefig6); 
        savefig(figure(7),namefig7);
        savefig(figure(8),namefig8);
        savefig(figure(9),namefig9);
    end 
end

% close all
% clear
end 



%% OTHER OPTION FOR HANDLING DATA FROM GIULIO
% 
% % After reshapping, row = mice, column = time bins
% ReshapedValues = transpose(reshape(in_Contact(:,1),[],[size(in_Contact,1)/numel(t_BinSizeIndices)]));
% ReshapedMouse = transpose(reshape(in_Contact(:,2),[],[size(in_Contact,1)/numel(t_BinSizeIndices)]));
% ReshapedTime = transpose(reshape(in_Contact(:,3),[],[size(in_Contact,1)/numel(t_BinSizeIndices)]));
% 
% [~,Order] = sort(ReshapedMouse(:,1),'ascend') ; % Reshape by pharmacological group 
% 
% ReshapedValues = ReshapedValues(Order,:);
% ReshapedMouse = ReshapedMouse(Order,:);
% ReshapedTime = ReshapedTime(Order,:);
% CumsumValues = cumsum(ReshapedValues,2) ;
% 
% figure(1);clf
%     plot(timeVector ,ReshapedValues(ReshapedMouse(:,1)==1,:),'r'); hold on,
%     plot(timeVector ,ReshapedValues(ReshapedMouse(:,1)==2,:),'k'); hold on,
%     plot(timeVector ,ReshapedValues(ReshapedMouse(:,1)==3,:),'color',[0.5,0.5,0.5]); hold on,
%     xlabel('Time(min)')
%     ylabel('duration of contact (s)')
%     title([num2str(Type_event),'- individual data'])
%     namefig1=[num2str(Type_event),'- individual data']; 
%     
% figure(2);clf    
%     plot(timeVector,cumsum(ReshapedValues(ReshapedMouse(:,1)==1,:),2),'r'); hold on ;
%     plot(timeVector,cumsum(ReshapedValues(ReshapedMouse(:,1)==2,:),2),'k'); hold on ;
%     plot(timeVector,cumsum(ReshapedValues(ReshapedMouse(:,1)==3,:),2),'color',[0.5,0.5,0.5]); hold on ;
%     xlabel('Time(min)')
%     ylabel('cumulative time spent in contact (s)')
%     title([num2str(Type_event),'- individual data']) 
%     namefig2=[num2str(Type_event),'- individual data cumu'];
%     
% figure(4);clf    
%     errorbar(timeVector,nanmean(cumsum(ReshapedValues(ReshapedMouse(:,1)==1,:),2),1),...
%         nanstd(cumsum(ReshapedValues(ReshapedMouse(:,1)==1,:),2),1)./sqrt(Group_N(1,1)-1),'r'); hold on
%     errorbar(timeVector,nanmean(cumsum(ReshapedValues(ReshapedMouse(:,1)==2,:),2),1),...
%         nanstd(cumsum(ReshapedValues(ReshapedMouse(:,1)==2,:),2),1)./sqrt(Group_N(1,2)-1),'k'); hold on
%     errorbar(timeVector,nanmean(cumsum(ReshapedValues(ReshapedMouse(:,1)==3,:),2),1),...
%         nanstd(cumsum(ReshapedValues(ReshapedMouse(:,1)==3,:),2),1)./sqrt(Group_N(1,3)-1),'color',[0.5,0.5,0.5]); hold on
%     xlabel('Time(min)')
%     ylabel('duration of contact (s)')
%     title([num2str(Type_event),'- Mean data +/- SEM'])    
%     namefig4=[num2str(Type_event),'- Mean data'];