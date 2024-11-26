clear all;close all;

cd D:\Nextcloud\Home\Behaviour\PatternDiscrimination\Bee_Pattern_vs_Colour\col_vs_pattern

%compare single and multiparametric learning cuves
filenames={'blue_orange_learningCurves.xlsx','orange_blue_learningCurves_shape.xlsx','orange_blue_learningCurves_col.xlsx',...
'blue_turquoise_learningCurves.xlsx','blue_turquoise_learningCurves_shape.xlsx','teal_blue_learningCurves_col.xlsx','pattern_learningCurves_col.xlsx','shape_learningCurves_col.xlsx'};

%compare all pattern combis
% filenames={'yellow_blue_learningCurves.xlsx','blue_orange_learningCurves.xlsx','yellow_orange_learningCurves.xlsx','blue_turquoise_learningCurves.xlsx'};
% type='pattern_'

%compare all shape combis
% filenames={'yellow_blue_learningCurves_shape.xlsx','orange_blue_learningCurves_shape.xlsx','orange_yellow_learningCurves_shape.xlsx','blue_turquoise_learningCurves_shape.xlsx'};
% type='shape_'

%compare all single attributes
% filenames={'orange_blue_learningCurves_col.xlsx','teal_blue_learningCurves_col.xlsx','pattern_learningCurves_col.xlsx','shape_learningCurves_col.xlsx'};
% type='single_'

if ~exist('filenames')
filenames=uigetfile('*learningCurves*.xlsx','multiselect','on');
end
%if only one file is selected
if iscell(filenames)==0
    temp=filenames; clear filenames;
    filenames=cell(1,1);
    filenames{1}=temp; clear temp;
end

allLastBlock=[];
groups=[];
allBlocks=cell(6,1); %collects learning data from all blocks

%initialise figure
f1=figure('Position',[500 500 700 250]);

for i=1:length(filenames)

    %% extract data

    %extract conditions from filename
    endInd=strfind(filenames{i},'_learning');
    colCond{i}=filenames{i}(1:endInd(1)-1);

    %extract first learning curve
    temp=readtable(filenames{i},"sheet",1,'ReadVariableNames',true);
    firstCurve{i}=temp{:,:};

    %extract all choice data from the last block
    allLastBlock=[allLastBlock;firstCurve{i}(6,:)'];
    groups=[groups;i*ones(size(firstCurve{i}(6,:)'))];

    for u=1:6
    allBlocks{u}=[allBlocks{u};firstCurve{i}(u,:)'];
    end

    %extract second learning curve
    temp=readtable(filenames{i},"sheet",2,'ReadVariableNames',true);
    secondCurve{i}=temp{:,:};

    %% plot figure
    subplot(1,2,1);hold on;
    mean=nanmean(firstCurve{i},2);sem=nanstd(firstCurve{i},1,2)/sqrt(size(firstCurve{i},2)-1);
%     boxplot(firstCurve{i}','PlotStyle','compact')
        errorbar(mean,sem,'.-','MarkerSize',10);
    xlabel('block')
    ylabel('fraction of correct choices')
    ylim([0 1])
    xlim([.5 6.5])

    subplot(1,2,2);hold on;
    mean=nanmean(secondCurve{i},2);sem=nanstd(secondCurve{i},1,2)/sqrt(size(secondCurve{i},2)-1);
        errorbar(mean,sem,'.-','MarkerSize',10);
%     boxplot(firstCurve{i}','PlotStyle','compact')

    xlabel('block')
    ylim([0 1])
    xlim([.5 6.5])
end

legend(colCond,"Location","south",'Interpreter','none')

%save figure
% print(f1,'learningCurves.eps','-dpdf','-r300','-painters','-bestfit')
% print(f1,'learningCurves_shape.eps','-dpdf','-r300','-painters','-bestfit')
% print(f1,'singleParams.eps','-dpdf','-r300','-painters','-bestfit')

%% compare learning success in the last blocks

 f2=figure("Position",[500 500 800 300]);
    subplot(1,3,1);hold on;
    boxplot(allLastBlock,groups);
%     xlim([0.5 4.5]) %force layout for 4 categories
    plot(get(gca,'XLim'),0.5*ones(length(get(gca,'XLim'))),'k--');hold on

    for u=1:max(groups)
    n=numel(allLastBlock(groups==u));
%     plot(spread*rand(n,1)+u*ones(n,1)-0.5*spread,colConflChoice(groups==u)/10,'.','MarkerSize',8,'color','k')
    s=swarmchart(u*ones(n,1),allLastBlock(groups==u),'.','k');
    s.XJitterWidth = 0.75;
%     b=beeswarm(u*ones(n,1),colConflChoice(groups==u)/10,'sort_style','hex','use_current_axes',true,'dot_size',0.2,'MarkerFaceColor','k','MarkerFaceAlpha',1);
    %analyse stats: difference to 0.5 (random choice)
    if sum(~isnan(allLastBlock(groups==u)))>0
    [p_confl(u),H,stats_confl] = signrank(allLastBlock(groups==u),5,"tail","right");
    signedrank_confl(u)=stats_confl.signedrank;
    else
        p_confl(u)=nan;signedrank_confl(u)=nan;
    end
    end

    %plot stats
    plot(get(gca,'XTick'),0.1*(p_confl<0.05),'*')

    set(gca,'XTickLabel',colCond)
    ylabel('fraction of correct choices')
    ylim([0 1])

    %save data for further statistical tests
    integerChoices=round(5./allLastBlock);%for learning curves, each block has always 5 correct choices, and x number of incorrect
    t=table([1:length(allLastBlock)]',5*ones(size(integerChoices)),integerChoices-5,groups,'VariableNames',{'AnimalID','correctChoices','incorrectChoices','group'});
 writetable(t,'allLastBlock.xlsx');

     %save data for statistical tests across all blocks
%      cd blockTests
%      for u=1:6
%          integerChoices=round(5./allBlocks{u});%for learning curves, each block has always 5 correct choices, and x number of incorrect
%     t=table([1:length(allBlocks{u})]',5*ones(size(integerChoices)),integerChoices-5,groups,'VariableNames',{'AnimalID','correctChoices','incorrectChoices','group'});
%  writetable(t,['allBlock_',type,num2str(u),'.xlsx']);
%      end
%     