% clear;
% close all;
% 

load('iML1515_lyco.mat')
% model=iML1515;

% get biomass reaction
model.csense=columnVector(model.csense);
biomassRxn=model.rxns{model.c==1};


% define target produc
targetRxn='EX_lyco_e';
targetID=findRxnIDs(model,targetRxn);

% set uptake rate of oxygen and glucose
oxygenRxn='EX_o2_e';
substrate='EX_glc__D_e';
model = changeRxnBounds(model,{substrate,oxygenRxn},-20,'l');

solWT=optimizeCbModel(model);

orimodel=model;

% when the model size is big, it is better to compress the model so that
% the linear reactions can be reduced
if size(model.S,1)>1000
    % compress the model and get compressed candidate reactions for knockout
    [model,candidate]=nihba_prep(orimodel,substrate,oxygenRxn,biomassRxn,targetRxn);
else
    candidate.rxns=model.rxns;
end

% limit reaction rate in realistic range
model.lb(model.lb<-100)=-100;
model.ub(model.ub>100)=100;

% make sure target reaction is in the reaction list of compressed model
if ~strcmp(model.rxns, targetRxn)
    targetRxn=model.rxns{contains(model.rxns, targetRxn)};
end

% remove unreasonable reactions from candidate knockout set
selectedRxns=setdiff(candidate.rxns, {'ATPM', biomassRxn, targetRxn});

% minimum ratio of growth in production mutants
minRationOfGrowth=0.0; % 10% wild type growth
model.lb(model.c==1)=minRationOfGrowth*solWT.f;

% call seqMOMAKnock and save results in excel and mat files
K=10;
mutant=seqMOMAKnock(model, selectedRxns,targetRxn, K);
save('moma-result-lyco.mat','mutant')