%% This code is an implementation of a sequential MOMA approach for identifying 
% multiple knockout solutions, which can be found in the following study:
% Alper H. et al. (2005), Identifying gene targets for the metabolic
% engineering of lycopene biosynthesis in Escherichia coli, Metabolic
% Engineering 7(2005) 155-164
function mutant=seqMOMAKnock(model, reactionList,targetRxn, K)
% INPUT:
%       model        - a cobra model structure
%       reactionList - a list of reations that can be considered as knockout targets
%       targetRxn    - biochemical of interest for overproduction
%       K            - the maximum number of knockouts allowed in a knockout strategy
% OUTPUT:
%       mutant       - identified knockout strategies with a structure that includes
%                    *.key   -storing knockout strategies for each number of knockouts
%                    *.value - storing growth-production pair for each strategy

% NOTE: You can use mutant(k).key (or mutant(k).value) to get designed strategies 
%(or growth-production pair) with k knockouts.

targetID=findRxnIDs(model,targetRxn);

solWT0=optimizeCbModel(model);

if ~iscell(reactionList)
    warning('>>Candidate knockout set is not cell contained.\n')
    reactionList={reactionList};
end

% param.Presolve=0;
param.OutputFlag=0;
param.OptimalityTol=1e-9;
param.FeasibilityTol=1e-9;
%
mutant=[];
candidateSet={};

%% section of single reaction knockout
k=1;
mutant(k).key={};
mutant(k).value={};
for i=1:length(reactionList)
    tmpRxnKO=reactionList{i};
    modelDel=changeRxnBounds(model,tmpRxnKO,0, 'b');
    try
        [solutionDel, solutionWT, totalFluxDiff] = MOMA(model, modelDel, 'max', 0,1);
        growth_target=[solutionDel.x(model.c==1), solutionDel.x(targetID)];
        
%         QCPproblem=buildQCP(modelDel, solutionWT.x, totalFluxDiff);
%         resultgurobi = gurobi(QCPproblem,param);
%         growth_target=[resultgurobi.x(model.c==1), resultgurobi.x(targetID)];
%         % keep nonzero production candidates
        if growth_target(2)>1e-8 && growth_target(1)>0.1*solWT0.f
            mutant(k).key=[mutant(k).key;tmpRxnKO];
            mutant(k).value=[mutant(k).value;growth_target];
            if k==1
                candidateSet{end+1}= tmpRxnKO;
            end
        end
        
%         if growth_target(2)>2.5
%                     a=1;
%         end
    catch
        warning('this mutant is not viable')
    end
end

%% section of multi-reaction knockout strategies
for k=2:K
   mutant(k).key={};
   mutant(k).value={}; 
   for i=1:length(mutant(k-1).key)
       koset=mutant(k-1).key(i,:);
       if ~iscell(koset) koset={koset}; end % single-ko
       modelWT=changeRxnBounds(model,koset, 0);
       
       [~,loc]=ismember(koset{end},candidateSet);
       for j=loc+1:length(candidateSet)
           tmpRxnKO=candidateSet{j};
           modelDel=changeRxnBounds(modelWT,tmpRxnKO, 0, 'b');
            try
                [solutionDel, solutionWT, totalFluxDiff] = MOMA(modelWT, modelDel, 'max', 0,true);
                growth_target=[solutionDel.x(model.c==1), solutionDel.x(targetID)];
                
%                 QCPproblem=buildQCP(modelDel, solutionWT.x, totalFluxDiff);
%                 resultgurobi = gurobi(QCPproblem,param);
%                 growth_target=[resultgurobi.x(model.c==1), resultgurobi.x(targetID)];
                
                
                % keep nonzero production candidates
                if growth_target(1)>0.1*solWT0.f && (growth_target(2)>1e-3 && growth_target(2)>=mutant(k-1).value{i}(2) && growth_target(2)>=mutant(1).value{j}(2))
                    mutant(k).key(end+1,:)=[koset,tmpRxnKO];
                    mutant(k).value{end+1}=growth_target;

                end
            catch
                warning('this mutant is not viable')
            end
       end
   end
   mutant(k).value=columnVector(mutant(k).value); % shown as column vectors
   save(['moma-result-',num2str(k),'.mat'],'mutant');
end

% % put all mutants together
mutants=[];
mutants.key={};
mutants.value={};
for k=1:K
    kk=length(mutant(k).key);
    mutants.key(end+1:end+kk)=mutant(k).key(1:kk);
    mutants.value(end+1:end+kk)=mutant(k).value(1:kk);
end
    
end


% CreateQCP that maximise growth subject to MOMA
function QCPproblem=buildQCP(model, w, delta)
    nRxns=length(w);
    QCPproblem = buildLPproblemFromModel(model);
    
    % Add second-order cone: x^2 + y^2 <= z^2
    QCPproblem.quadcon(1).Qc = speye(nRxns);
    QCPproblem.quadcon(1).q=-2*w';
    QCPproblem.quadcon(1).rhs=delta+min(1e-9, 1e-9*delta)-dot(w',w);
    QCPproblem.quadcon(1).name = 'moma';
end