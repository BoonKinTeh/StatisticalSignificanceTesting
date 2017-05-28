function Result = SignificanceTesting_PL(Data,PerformInvFormale,PerformCSN)
%% Inputs description:
% Data: Is a 1xN array contains N observations (hypothetically power law distributed)
% PerformInvFormale: set to 1 if wish to perform the fast inversion formule of statistical significance testing
% PerformCSN: set to 1 isf wish to perform the CSN algorithm of statistical significance testing (Very slow process, O(size(Data)*CSN_NSim*XminANumPartition))
%% Default inputs description:
CSN_NSim = 500; % Number of bootstrapped samples for statistical significance testing
XminANumPartition = 500; % Length of the Xmin Array (For speed up purposes), set to max(size(Data)) if don't wish to speed up.
                         % Change the power law spaced Xmin array to linear spaced Xmin array
%% Outputs description:
% Result: A structure format data, contain information about estimated parameters, 
%         Inversion Formale statistical significance testing (if PerformInvFormale==1),
%         and CSN statistical significance testing (if PerformCSN ==1)
% To call output: Output = Result.Output
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
%% Initialization 
Data = reshape(Data,[1,max(size(Data))]); % Make sure Data is 1xN array
Data = sort(Data);

Sample_Alpha = 0;
Sample_Xmin = 0;
Sample_KSDistance = 0;
Sample_DistributionNoise = 0;
Sample_NumFittedSample = 0;

XminA = unique(Data); XminA = reshape(XminA,[1,size(XminA,1)*size(XminA,2)]);
Ind = XminA>=min(Data(:,end+1-min(size(Data,2),50)));  XminA(:,Ind) = [];
if size(XminA,2)>XminANumPartition
    XminA = linspace(min(XminA),max(XminA),XminANumPartition);
end
PreDist = inf;

%% Parameter Estimations
warning('off','all');
options = optimoptions('lsqnonlin','Display','off');
for Xmin_i = 1:size(XminA,2)-1
    Xmin = XminA(1,Xmin_i);
    Sample = Data(Data>=Xmin);
    NumFittedSample = size(Sample,2);
    Avg_LogXPerXmin = mean(log(Sample))-log(Xmin);
    %%% Estimate Alpha
    Alpha = 1+1./Avg_LogXPerXmin; XMax = max(Sample);
    AlphaFunction =@(A) abs( (Avg_LogXPerXmin+1./(1-A))+ (log(XMax./Xmin)+1./(A-1)-Avg_LogXPerXmin).*(XMax./Xmin).^(1-A) );
    Alpha = lsqnonlin(AlphaFunction,Alpha,2,inf,options);
    %%% Probability Integral Transform
    CDF_Th = linspace(0,1,NumFittedSample);
    Delta = (XMax./Xmin).^(1-Alpha);
    CDF_Empirical = 1-(Sample./Xmin).^(1-Alpha);
    CDF_Empirical = CDF_Empirical./(1-Delta);
    %%% KS Distance Measure
    KSDistance = max(abs(CDF_Th-CDF_Empirical));    
    %%% Distribution Noise Measure
    Spacing = diff([0,CDF_Empirical]);  Spacing(:,Spacing==0)=[];
    DistributionNoise = sqrt(sum((Spacing.^2).*(1/NumFittedSample./Spacing-1).^2)/sum(Spacing.^2));
    %%% Estimate Xmin (Shortest KS Distance)
    Sample_Alpha = (PreDist<=KSDistance)*Sample_Alpha + (PreDist>KSDistance)*Alpha;
    Sample_Xmin = (PreDist<=KSDistance)*Sample_Xmin + (PreDist>KSDistance)*Xmin;
    Sample_KSDistance = (PreDist<=KSDistance)*Sample_KSDistance + (PreDist>KSDistance)*KSDistance;
    Sample_DistributionNoise = (PreDist<=KSDistance)*Sample_DistributionNoise + (PreDist>KSDistance)*DistributionNoise;
    Sample_NumFittedSample = (PreDist<=KSDistance)*Sample_NumFittedSample + (PreDist>KSDistance)*NumFittedSample;
    PreDist = min(PreDist,KSDistance);
end
Result.Data = Data;
Result.Sample_Alpha = Sample_Alpha;
Result.Sample_Xmin = Sample_Xmin;
Result.Sample_KSDistance = Sample_KSDistance;
Result.Sample_DistributionNoise = Sample_DistributionNoise;
Result.Sample_NumFittedSample = Sample_NumFittedSample;
clearvars -except Data PerformInvFormale PerformCSN XminANumPartition CSN_NSim Result...
            Sample_Alpha Sample_Xmin Sample_KSDistance Sample_DistributionNoise Sample_NumFittedSample 

%% Statistical Significance Testing: Fast Inversion Formule (FIF)
if PerformInvFormale==1
    options = optimoptions('lsqnonlin','Display','off','TolFun',1e-12);
    ParaKS = [0.2744,0.1760,0.4924];
    ParaDN = [0.4298,0.3023,0.4946];
    FIF_KS = 100./(1+(exp(ParaKS(1,1))*Sample_KSDistance.* (Sample_NumFittedSample.^(ParaKS(1,3)))).^(-1/ParaKS(1,2)));

    EDN = (Sample_NumFittedSample/(0.5+Sample_NumFittedSample))*sqrt(1/2+ (2-Sample_NumFittedSample)/(2*Sample_NumFittedSample^2) );
    Miu = Sample_DistributionNoise-EDN;
    if Miu<0
        PDNFun =@(PDN)   PDN.^(ParaDN(1,1))+ log(abs(Miu).*(Sample_NumFittedSample.^ParaDN(1,3))).*((50-PDN).^(ParaDN(1,2))) ;
        FIF_DN = lsqnonlin(PDNFun,25,0,50,options);
    elseif Miu==0
        FIF_DN = 50;
    else
        PDNFun =@(PDN) (100-PDN).^(ParaDN(1,1))+ log(abs(Miu).*(Sample_NumFittedSample.^ParaDN(1,3))).*((PDN-50).^(ParaDN(1,2)));
        FIF_DN = lsqnonlin(PDNFun,75,50,100,options);
    end
    FIF_KSDN = 100*sqrt((1-FIF_KS/100)*(1-FIF_DN/100)*(1-exp(1)/(Sample_NumFittedSample^0.4809)));
    FIF_KS = 100-FIF_KS;
    FIF_DN = 100-FIF_DN;
    
    Result.FIF_KS = FIF_KS;
    Result.FIF_DN = FIF_DN;
    Result.FIF_KSDN = FIF_KSDN;
    clearvars -except Data PerformInvFormale PerformCSN XminANumPartition CSN_NSim Result...
            Sample_Alpha Sample_Xmin Sample_KSDistance Sample_DistributionNoise Sample_NumFittedSample 
end

%% Statistical Significance Testing: Fast Inversion Formule (FIF)
if PerformCSN==1
    CSN__KSDistance = zeros(CSN_NSim,1);
    CSN__DistributionNoise = zeros(CSN_NSim,1);
    CSN__NumFittedSample = zeros(CSN_NSim,1);
    rng('shuffle');
    warning('off','all');
    options = optimoptions('lsqnonlin','Display','off');
    for Run_i = 1:CSN_NSim
        %%% Bootstrapping for CSN Significance Testing
        SegmentBeforeXmin = Data(Data<Sample_Xmin);
        SegmentBeforeXmin = SegmentBeforeXmin(1,ceil(size(SegmentBeforeXmin,2)*rand(1,size(SegmentBeforeXmin,2))));
        SampleSize = size(Data,2)-size(SegmentBeforeXmin,2);
        ResampledSample  = Sample_Xmin*(rand(1,SampleSize)).^(1/(1-Sample_Alpha));
        BootstrappedSample = sort([SegmentBeforeXmin,ResampledSample]);
        %%% Repeat same process as sample parameter estimation
        XminA = unique(BootstrappedSample); XminA = reshape(XminA,[1,size(XminA,1)*size(XminA,2)]);
        Ind = XminA>=min(BootstrappedSample(:,end+1-min(size(BootstrappedSample,2),50)));  XminA(:,Ind) = [];
        if size(XminA,2)>XminANumPartition
            XminA = linspace(min(XminA),max(XminA),XminANumPartition);    
        end
        PreDist = inf;
        for Xmin_i = 1:size(XminA,2)-1
            Xmin = XminA(1,Xmin_i);
            Sample = BootstrappedSample(BootstrappedSample>=Xmin);
            NumFittedSample = size(Sample,2);
            Avg_LogXPerXmin = mean(log(Sample))-log(Xmin);
            %%% Estimate Alpha
            Alpha = 1+1./Avg_LogXPerXmin; XMax = max(Sample);
            AlphaFunction =@(A) abs( (Avg_LogXPerXmin+1./(1-A))+ (log(XMax./Xmin)+1./(A-1)-Avg_LogXPerXmin).*(XMax./Xmin).^(1-A) );
            Alpha = lsqnonlin(AlphaFunction,Alpha,2,inf,options);
            %%% Probability Integral Transform
            CDF_Th = linspace(0,1,NumFittedSample);
            Delta = (XMax./Xmin).^(1-Alpha);
            CDF_Empirical = 1-(Sample./Xmin).^(1-Alpha);
            CDF_Empirical = CDF_Empirical./(1-Delta);
            %%% KS Distance Measure
            KSDistance = max(abs(CDF_Th-CDF_Empirical));    
            %%% Distribution Noise Measure
            Spacing = diff([0,CDF_Empirical]);  Spacing(:,Spacing==0)=[];
            DistributionNoise = sqrt(sum((Spacing.^2).*(1/NumFittedSample./Spacing-1).^2)/sum(Spacing.^2));
            %%% Estimate Xmin (Shortest KS Distance)
            CSN__KSDistance(Run_i,1) = (PreDist<=KSDistance)*CSN__KSDistance(Run_i,1) + (PreDist>KSDistance)*KSDistance;
            CSN__DistributionNoise(Run_i,1) = (PreDist<=KSDistance)*CSN__DistributionNoise(Run_i,1) + (PreDist>KSDistance)*DistributionNoise;
            CSN__NumFittedSample(Run_i,1) = (PreDist<=KSDistance)*CSN__NumFittedSample(Run_i,1) + (PreDist>KSDistance)*NumFittedSample;
            PreDist = min(PreDist,KSDistance);
        end
    end
    KS_NumFittedSamplePower = 0.4924;
    DN_NumFittedSamplePower = 0.4946;
    
    CSN_KS = 100-100*sum(CSN__KSDistance.*(CSN__NumFittedSample.^KS_NumFittedSamplePower)<Sample_KSDistance.*(Sample_NumFittedSample.^KS_NumFittedSamplePower))/CSN_NSim;
    CSN_DN = 100-100*sum((CSN__DistributionNoise-1/sqrt(2)).*(CSN__NumFittedSample.^DN_NumFittedSamplePower)<(Sample_DistributionNoise-1/sqrt(2)).*(Sample_NumFittedSample.^DN_NumFittedSamplePower))/CSN_NSim;
    CSN_KSDN = 100*sqrt(CSN_KS/100*CSN_DN/100*(1-exp(1)/(Sample_NumFittedSample^0.4809)));
    
    Result.CSN_KS = CSN_KS;
    Result.CSN_DN = CSN_DN;
    Result.CSN_KSDN = CSN_KSDN;
end
