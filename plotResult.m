% this code is used plot the production envelope for mutant and wide-type in one figure
biomassRxn=model.rxns{model.c==1};
targetRxn='EX_succ_e';

% the following is to get the best strategy of k=5 knockouts with the highest production rate 
% from the sequential MOMA, and then visualise the production envelope
K=5;
allSet=mutant(K).key;
gr_prod=cell2mat(mutant(K).value);
[s_gp, idx]=sort(gr_prod(:,2),'descend');
best=allSet(idx(1),:);

% alternatively, you can manually specify a strategy you want forvisualisation
% best={'PGI','ATPS4rpp','PDH','EX_o2_e','LDH_D'};


model.csense=columnVector(model.csense);
geneDelFlag = false;
nPts = 30;
figure(1)
productionEnvelope(model,{},'r',targetRxn,biomassRxn,geneDelFlag,nPts);
xlabel('Biomass', 'FontSize', 20);
ylabel('Production Rate', 'FontSize', 20);
hold on
productionEnvelope(model,best, 'b',targetRxn,biomassRxn,geneDelFlag,nPts);