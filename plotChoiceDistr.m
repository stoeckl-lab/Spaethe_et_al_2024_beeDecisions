clear all

cd D:\Nextcloud\Home\Behaviour\PatternDiscrimination\Bee_Pattern_vs_Colour\col_vs_pattern

filenames=uigetfile('*Choices*.xlsx');

%check number of sheets (i.e. number of conditions)
[status,sheets] = xlsfinfo(filenames);
colCond=sheets;

%initialise variable
colConflChoice=[];
nonPrefPatChoice=[];
nonPrefColChoice=[];
groups=[];

for i=1:length(sheets)

    %% extract data

    %extract data from each condition learning curve
    temp=readtable(filenames,"sheet",i);
    if isempty(temp)
        temp=num2cell(nan(1,5));
    end

    colConflChoice=[colConflChoice;temp{:,2}];
    if size(temp,2)>2
    nonPrefPatChoice=[nonPrefPatChoice;temp{:,4}];
    nonPrefColChoice=[nonPrefColChoice;temp{:,5}];
    else
    nonPrefPatChoice=nan(size(colConflChoice));
    nonPrefColChoice=nan(size(colConflChoice));
    end
       groups=[groups;i*ones(size(temp{:,2},1),1)];
 
end

   %% plot figure
   spread=0.3;

    f1=figure("Position",[500 500 800 300]);
    subplot(1,3,1);hold on;
    boxplot(colConflChoice/10,groups);
%     xlim([0.5 4.5]) %force layout for 4 categories
    plot(get(gca,'XLim'),0.5*ones(length(get(gca,'XLim'))),'k--');hold on

    for u=1:max(groups)
    n=numel(colConflChoice(groups==u));
%     plot(spread*rand(n,1)+u*ones(n,1)-0.5*spread,colConflChoice(groups==u)/10,'.','MarkerSize',8,'color','k')
    s=swarmchart(u*ones(n,1),colConflChoice(groups==u)/10,'.','k');
    s.XJitterWidth = 0.75;
%     b=beeswarm(u*ones(n,1),colConflChoice(groups==u)/10,'sort_style','hex','use_current_axes',true,'dot_size',0.2,'MarkerFaceColor','k','MarkerFaceAlpha',1);
    %analyse stats: difference to 0.5 (random choice)
    if sum(~isnan(colConflChoice(groups==u)))>0
    [p_confl(u),H,stats_confl] = signrank(colConflChoice(groups==u),5,"tail","right");
    signedrank_confl(u)=stats_confl.signedrank;
    else
        p_confl(u)=nan;signedrank_confl(u)=nan;
    end
    end

    %plot stats
    plot(get(gca,'XTick'),0.1*(p_confl<0.05),'*')

    set(gca,'XTickLabel',colCond)
    ylabel('fraction of colour choices')
    ylim([0 1])
    

    subplot(1,3,2);hold on;
    boxplot(nonPrefPatChoice/10,groups);
    plot(get(gca,'XLim'),0.5*ones(length(get(gca,'XLim'))),'k--');hold on
    for u=1:max(groups)
    n=numel(nonPrefPatChoice(groups==u));
    s=swarmchart(u*ones(n,1),nonPrefPatChoice(groups==u)/10,'.','k');
    s.XJitterWidth = 0.75;
    %analyse stats: difference to 0.5 (random choice)
    if sum(~isnan(nonPrefPatChoice(groups==u)))>0
    [p_noprefPat(u),H,stats_noprefPat] = signrank(nonPrefPatChoice(groups==u),5,"tail","right");
    signedrank_noprefPat(u)=stats_noprefPat.signedrank;
    else
        p_noprefPat(u)=nan;signedrank_noprefPat(u)=nan;
    end

    end


     %plot stats
    plot(get(gca,'XTick'),0.1*(p_noprefPat<0.05),'*')

%     set(gca,'XTickLabel',colCond)
    ylabel('fraction of non pref pattern choices')
    ylim([0 1])

    subplot(1,3,3);hold on;
    boxplot(nonPrefColChoice/10,groups);
    plot(get(gca,'XLim'),0.5*ones(length(get(gca,'XLim'))),'k--');hold on
    for u=1:max(groups)
    n=numel(nonPrefColChoice(groups==u));
    s=swarmchart(u*ones(n,1),nonPrefColChoice(groups==u)/10,'.','k');
    s.XJitterWidth = 0.75;
    if sum(~isnan(nonPrefColChoice(groups==u)))>0
    %analyse stats: difference to 0.5 (random choice)
    
    [p_noprefCol(u),H,stats_noprefCol] = signrank(nonPrefColChoice(groups==u),5,"tail","right");
    signedrank_noprefCol(u)=stats_noprefCol.signedrank;
    else 
        p_noprefCol(u)=nan;signedrank_noprefCol(u)=nan;
    end
    end

     %plot stats
    plot(get(gca,'XTick'),0.1*(p_noprefCol<0.05),'*')
%     set(gca,'XTickLabel',colCond)

    ylabel('fraction of non pref col choices')
    ylim([0 1])

    %% save figure
% print(f1,'choiceDistribution.eps','-dpdf','-r300','-painters','-bestfit')
% print(f1,'choiceDistribution_shape.eps','-dpdf','-r300','-painters','-bestfit')

