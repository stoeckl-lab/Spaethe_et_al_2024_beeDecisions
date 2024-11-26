
filenames=uigetfile('*Choices.xlsx');

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
    %extract data from each condition 
    temp=readtable(filenames,"sheet",i);
    %always two columns belong together as the choices that were presented
    data=temp{:,2:5}/10;
    
    groups=[groups;i*ones(size(temp{:,2},1),1)];
 
end

   %% plot figure
   spread=0.2;

    figure("Position",[500 500 800 400]);
    subplot(1,3,1);hold on;
    boxplot(colConflChoice/10,groups);
    plot(get(gca,'XLim'),0.5*ones(length(get(gca,'XLim'))),'k--');hold on
    for u=1:max(groups)
    n=numel(colConflChoice(groups==u));
    plot(spread*rand(n,1)+u*ones(n,1)-0.5*spread,colConflChoice(groups==u)/10,'.','MarkerSize',12,'color','k')
    %analyse stats: difference to 0.5 (random choice)
    [p_confl(u),H,stats_confl] = signrank(colConflChoice(groups==u),5,"tail","right");
    signedrank_confl(u)=stats_confl.signedrank;
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
    plot(spread*rand(n,1)+u*ones(n,1)-0.5*spread,nonPrefPatChoice(groups==u)/10,'.','MarkerSize',12,'color','k')
    %analyse stats: difference to 0.5 (random choice)
    [p_noprefPat(u),H,stats_noprefPat] = signrank(nonPrefPatChoice(groups==u),5,"tail","right");
    signedrank_noprefPat(u)=stats_noprefPat.signedrank;
    end

     %plot stats
    plot(get(gca,'XTick'),0.1*(p_noprefPat<0.05),'*')

    set(gca,'XTickLabel',colCond)
    ylabel('fraction of non pref pattern choices')
    ylim([0 1])

    subplot(1,3,3);hold on;
    boxplot(nonPrefColChoice/10,groups);
    plot(get(gca,'XLim'),0.5*ones(length(get(gca,'XLim'))),'k--');hold on
    for u=1:max(groups)
    n=numel(nonPrefColChoice(groups==u));
    plot(spread*rand(n,1)+u*ones(n,1)-0.5*spread,nonPrefColChoice(groups==u)/10,'.','MarkerSize',12,'color','k')
    if sum(~isnan(nonPrefColChoice(groups==u)))>0
    %analyse stats: difference to 0.5 (random choice)
    [p_noprefCol(u),H,stats_noprefCol] = signrank(nonPrefColChoice(groups==u),5,"tail","right");
    signedrank_noprefCol(u)=stats_noprefCol.signedrank;
    else 
        p_noprefCol(u)=nan;
    end
    end
    
     %plot stats
    plot(get(gca,'XTick'),0.1*(p_noprefCol<0.05),'*')
set(gca,'XTickLabel',colCond)

    ylabel('fraction of non pref col choices')
    ylim([0 1])

    %% stats

