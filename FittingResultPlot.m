function FittingResultPlot(Result)
%% Inputs description:
% Result: A structure data contains the significance testing result.
%% Default inputs description:
XLabelName = 'x';
YLabel = 'CDF';
LabelFontsize = 14;     % Label font size
AxisFontsize = 14;          % Axis font size 
TextFontsize = 16;      % Font size for the statistical significance value
%% Read Me:
% This code is Published together with "Finite sample corrections for parameters estimation and significance testing", 
% by Boon Kin Teh, Darrell JiaJie Tay, Sai Ping Li, and Siew Ann Cheong.
% Please refer to the paper for more details, and cite the paper when you
% are using this code for significance testing analysis,
% Thank you.

%% Lastest updated date:
% 27 May 2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract Information
DistributionType = 1*isfield(Result, 'Sample_Alpha')+2*isfield(Result, 'Sample_Beta');
%%% Sample details, hence constructs empirical and fits 1-cummulative density function
X = sort(Result.Data);
Xmin = Result.Sample_Xmin;
X_Fits = linspace(Xmin,max(X),50);
if DistributionType==1 
    Alpha = Result.Sample_Alpha;
    delta = (max(X)./Xmin).^(1-Alpha);
    Delta = delta*sum(X>=Xmin)/(size(X,2)-delta*sum(X<Xmin));
    ICDF_Empirical = linspace(1,Delta,size(X,2));
    ICDF_Fits = ICDF_Empirical(sum(X<Xmin))*(X_Fits./Xmin).^(1-Alpha);
elseif DistributionType==2
    Beta = Result.Sample_Beta;
    delta = exp(-Beta.*(max(X)-Xmin));
    Delta = delta*sum(X>=Xmin)/(size(X,2)-delta*sum(X<Xmin));
    ICDF_Empirical = linspace(1,Delta,size(X,2));
    ICDF_Fits = ICDF_Empirical(sum(X<Xmin))*exp(-Beta*(X_Fits-Xmin));
end
%%% Significance testing details.
SignificanceTestingType = 1*isfield(Result, 'FIF_KS')+2*isfield(Result, 'CSN_KS');
if isfield(Result, 'FIF_KS')
    FIF_KS = Result.FIF_KS;
    FIF_DN = Result.FIF_DN;
    FIF_KSDN = Result.FIF_KSDN;
end
if isfield(Result, 'CSN_KS')
    CSN_KS = Result.CSN_KS;
    CSN_DN = Result.CSN_DN;
    CSN_KSDN = Result.CSN_KSDN;
end
if SignificanceTestingType==1
    Text_KS = ['P^{FIF}_{KS}    : ' sprintf('%0.1f',FIF_KS)];
    Text_DN = ['P^{FIF}_{DN}    : ' sprintf('%0.1f',FIF_DN)];
    Text_KSDN = ['P^{FIF}_{KSDN} : ' sprintf('%0.1f',FIF_KSDN)];
elseif SignificanceTestingType==2
    Text_KS = ['P^{CSN}_{KS}    : ' sprintf('%0.1f',CSN_KS)];
    Text_DN = ['P^{CSN}_{DN}    : ' sprintf('%0.1f',CSN_DN)];
    Text_KSDN = ['P^{CSN}_{KSDN} : ' sprintf('%0.1f',CSN_KSDN)];
elseif SignificanceTestingType==3
    Text_KS = ['P^{CSN}_{KS} / P^{FIF}_{KS}      : ' sprintf('%0.1f',CSN_KS) ' / ' sprintf('%0.1f',FIF_KS)];
    Text_DN = ['P^{CSN}_{DN} / P^{FIF}_{DN}      : ' sprintf('%0.1f',CSN_DN) ' / ' sprintf('%0.1f',FIF_DN)];
    Text_KSDN = ['P^{CSN}_{KSDN} / P^{FIF}_{KSDN} : ' sprintf('%0.1f',CSN_KSDN) ' / ' sprintf('%0.1f',FIF_KSDN)];
end
%% Prepare for Plot
%%% Xtick and Ytick and text postion
if DistributionType==1
    XTick = 10.^(floor(log(min(X))/log(10)):ceil(log(max(X))/log(10)));
    TextPostionX = 10.^linspace(log(XTick(1,1))/log(10),log(XTick(1,end))/log(10),250);
    XScale = 'log';
elseif DistributionType==2
    Digi = floor(log(max(X))/log(10));
    XTick =(floor(min(X)/(10^Digi)):ceil(max(X)/(10^Digi))).*10^Digi;
    TextPostionX = linspace(XTick(1,1),XTick(1,end),250);
    XScale = 'linear';
end
YTick = 10.^(floor((log(ICDF_Empirical(1,end))/log(10))):ceil((log(1)/log(10))));
TextPostionY = 10.^linspace((log(YTick(1,1))/log(10)),(log(YTick(1,end))/log(10)),250);
%% Plotting
figure(1);clf;hold on;
plot(X,ICDF_Empirical,'.k','markersize',10);
plot(X_Fits,ICDF_Fits,'--','color',[1 0 0],'linewidth',3);
Legend = legend('Empirical','Fits','Location','bestoutside','Orientation','horizontal');
set(Legend,'FontSize',AxisFontsize);
plot(X(sum(X<Xmin)),ICDF_Fits(1,1),'h','color',[1 0 0],'MarkerFaceColor',[1 0 0],'markersize',15);
xlim([XTick(1,1),XTick(1,end)]); 
ylim([YTick(1,1),YTick(1,end)]);
set(gca,'xtick',XTick,'ytick',YTick);
set(gca,'xscale',XScale,'yscale','log','fontsize',AxisFontsize);
if SignificanceTestingType>0
    text(TextPostionX(1,20),TextPostionY(1,75),Text_KS,'color','r','fontsize',TextFontsize);
    text(TextPostionX(1,20),TextPostionY(1,50),Text_DN,'color','b','fontsize',TextFontsize);
    text(TextPostionX(1,20),TextPostionY(1,25),Text_KSDN,'color',[0,0.75,0.25],'fontsize',TextFontsize);
    text(TextPostionX(1,end-50),TextPostionY(1,end-10),['\delta: ' sprintf('%0.3f',delta)],'color','b','fontsize',2+TextFontsize);
end
xlabel(XLabelName,'fontsize',LabelFontsize); ylabel(YLabel,'fontsize',LabelFontsize);

