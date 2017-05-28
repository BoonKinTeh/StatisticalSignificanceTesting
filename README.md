# Statistical Significance Testing
This project is published with "Finite sample corrections for parameters estimation and significance testing", 
by Boon Kin Teh, Darrell JiaJie Tay, Sai Ping Li, and Siew Ann Cheong. Please refer to the paper for more details, and cite the paper if you are using this code for significance testing analysis, Thank you.

###########################################################################################
Code: SignificanceTesting_PL
Perform statistical significance testing for power law distribution. 
Inputs: 
1) Data: 1xN array for observations (hypothetically power law distributed)
2) PerformInvFormale: Set to 1 if wish to perform the significance testing using fast inversion formule (Fast).
3) PerformCSN: Set to 1 if wish to perform the significance testing using CSN method (very slow).
Outputs:
Result = SignificanceTesting_PL(Data,1,1);
1) Result: A structure format data, to call output: 
Data = Result.Data;   %Inputs (1xN array for observations)
Sample_Alpha = Result.Sample_Alpha;   %Estimated Alpha
Sample_Xmin = Result.Sample_Xmin;   %Estimated Xmin
Sample_KSDistance = Result.Sample_KSDistance;     %Measured KS distance
Sample_DistributionNoise = Result.Sample_DistributionNoise;     %Measured Distribution noise
Sample_NumFittedSample = Result.Sample_NumFittedSample;     %Number of data fitted to power law distribution
%(if PerformInvFormale==1)
FIF_KS = Result.FIF_KS;     %P-value (in percentage) based on KS distance using Fast inversion formule
FIF_DN = Result.FIF_DN;     %P-value (in percentage) based on distribution noise using Fast inversion formule
FIF_KSDN = Result.FIF_KSDN;     %P-value (in percentage) based on KS distance and distribution noise using Fast inversion formule
%(if PerformCSN==1)
CSN_KS = Result.CSN_KS;     %P-value (in percentage) based on KS distance using CSN method
CSN_DN = Result.CSN_DN;     %P-value (in percentage) based on distribution noise using CSN method
CSN_KSDN = Result.CSN_KSDN;     %P-value (in percentage) based on KS distance and distribution noise using CSN method


###########################################################################################
Code: SignificanceTesting_EXP
Perform statistical significance testing for exponential distribution. 
Inputs: 
1) Data: 1xN array for observations (hypothetically exponential distributed)
2) PerformInvFormale: Set to 1 if wish to perform the significance testing using fast inversion formule (Fast).
3) PerformCSN: Set to 1 if wish to perform the significance testing using CSN method (very slow).
Outputs:
Result = SignificanceTesting_PL(Data,1,1);
1) Result: A structure format data, to call output: 
Data = Result.Data;   %Inputs (1xN array for observations)
Sample_Alpha = Result.Sample_Alpha;   %Estimated Alpha
Sample_Xmin = Result.Sample_Xmin;   %Estimated Xmin
Sample_KSDistance = Result.Sample_KSDistance;     %Measured KS distance
Sample_DistributionNoise = Result.Sample_DistributionNoise;     %Measured Distribution noise
Sample_NumFittedSample = Result.Sample_NumFittedSample;     %Number of data fitted to power law distribution
%(if PerformInvFormale==1)
FIF_KS = Result.FIF_KS;     %P-value (in percentage) based on KS distance using Fast inversion formule
FIF_DN = Result.FIF_DN;     %P-value (in percentage) based on distribution noise using Fast inversion formule
FIF_KSDN = Result.FIF_KSDN;     %P-value (in percentage) based on KS distance and distribution noise using Fast inversion formule
%(if PerformCSN==1)
CSN_KS = Result.CSN_KS;     %P-value (in percentage) based on KS distance using CSN method
CSN_DN = Result.CSN_DN;     %P-value (in percentage) based on distribution noise using CSN method
CSN_KSDN = Result.CSN_KSDN;     %P-value (in percentage) based on KS distance and distribution noise using CSN method


###########################################################################################
Code: FittingResultPlot
Plot the fitting results
Inputs: 
1) Result: Result obtained from SignificanceTesting_PL or SignificanceTesting_EXP
Outputs:
A plot for the fitting results
