function c=config()
% Fixed study settings; all paths are relative to this package.
c.root=fileparts(mfilename('fullpath'));
addpath(fullfile(c.root,'algorithms'),fullfile(c.root,'functions'));
c.networks={'Ctown','Modena'};c.conditions={{'H1','H2'},{'M1','M2'}};
c.runs=20;c.population=30;c.FE=struct('Ctown',4000,'Modena',4000);
c.main={'EPPE','PPE','EROA','SSA','SOA','DBO'};
c.ablation={'PPE','PPE_M1','PPE_M2','EPPE','RandomExchange'};
c.seedBase=2520110101;c.statsSeed=2530110101;
c.stagnation=4;c.recipientFraction=.2;c.twoSwapProbability=.5;c.exploreProbability=.2;
c.B=struct('Ctown',21,'Modena',14);c.existing=struct('Ctown',4,'Modena',3);
c.floorSlack=.02;c.feasibilityTolerance=1e-10;
c.forbiddenFraction=.2;c.priorityFraction=.4;c.partitionSeedBase=740010;c.maskSeedBase=2510110101;
c.threshold=struct('Ctown',[2000 20 .8],'Modena',[2000 20 .8]);
c.accuracy=1e-5;c.maxTrials=200;c.pressureTolerance=1e-6;c.maxInvalidFraction=.05;c.impulseFraction=.1;
c.normalMILPseconds=60;c.dataDir=fullfile(c.root,'results','hydraulics');
c.toolkit='';
end
