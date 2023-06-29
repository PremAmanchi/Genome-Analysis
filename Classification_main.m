%classification using total length of sequences, 10 fold, using complexities(ETC/LZ/ETC+LZ) 
%as features for four classifiers(linear SVM(LSVM), Quadratic SVM(QSVM), Linear Discriminant(LD), Fine KNN(FKNN)) 
clc
clear all
close all
warning off

dataSet = 'covid';
selectedFolder = dataSet;
fprintf('Reading sequences .... \n');
[AcNmb, Seq, numberOfClusters, clusterNames, pointsPerCluster] = readFasta(dataSet);
totalSeq = length(Seq);

%calculate length stats
[maxLen, minLen, meanLen, medLen] = lengthCalc(Seq);

Y=[];
for i=1:numberOfClusters
    for j=1:pointsPerCluster{i}
        Y=[Y; i];
    end
end
classes=unique(Y);
ms=length(classes);

feature=input("Enter complexity measure(ETC/LZ/ETC+LZ): ");

if (feature=="ETC")
%{
ETC complexity calculation
nmValSH=cell(1,totalSeq);
lg=cell(1,totalSeq);
noretc=cell(1,totalSeq);
fprintf('Generating numerical sequences, applying ETC .... \n');
parfor a = 1:totalSeq  %parallel for loop
    ns = numMappingPP(Seq{a});
    %change "medLen" to other length stat for length normalization
    nmValSH{a} = ns;
    lg{a} = ETC(ns,0);
    noretc{a}=lg{a}/(length(Seq{a})-1);
end
X=cell2mat(noretc);
X=X';
%}
    
X=readmatrix('covid_etc_values.xlsx');

elseif (feature=="LZ")
%{
LZ complexity calculation
lg=cell(1,totalSeq);
fprintf('applying Lempel-ziv for sequences .... \n');
for a = 1:totalSeq  %parallel for loop
    
    [codice, code_bin, code_book]= lempel_ziv(['A' 'T' 'C' 'G'],cell2mat(Seq(a)));
    lg{a}=length(code_book);
     
    norlem{a}=(lg{a}/(length(Seq{a})))*(log(length(Seq{a}))/log(4));
       
end
X=cell2mat(norlem);
X=X';
%}
X=readmatrix('covid_lz_values.xlsx');

elseif (feature=="ETC+LZ") 
%ETC+LZ complexity calculation
nmValSH=cell(1,totalSeq);
lg=cell(1,totalSeq);
noretc=cell(1,totalSeq);
fprintf('Generating numerical sequences, applying ETC .... \n');
parfor a = 1:totalSeq  %parallel for loop
    ns = numMappingPP(Seq{a});
    %change "medLen" to other length stat for length normalization
    nmValSH{a} = ns;
    lg{a} = ETC(ns,0);
    noretc{a}=lg{a}/(length(Seq{a})-1);
end
lg=cell(1,totalSeq);
fprintf('applying Lempel-ziv for sequences .... \n');
for a = 1:totalSeq  %parallel for loop
    
    [codice, code_bin, code_book]= lempel_ziv(['A' 'T' 'C' 'G'],cell2mat(Seq(a)));
    lg{a}=length(code_book);
     
    norlem{a}=(lg{a}/(length(Seq{a})))*(log(length(Seq{a}))/log(4));
       
end
X=cell2mat(noretc);
X=X';
X1=cell2mat(norlem);
X1=X1';
X(:,2)=X1(:,1);

%X=readmatrix('covid_etc_values.xlsx');
%X1=readmatrix('covid_lz_values.xlsx');
%X(:,2)=X1(:,1);

else
    fprintf('Wrong feature provided! \n');
end

fprintf('Performing classification.... \n');
folds=10;
%linear discriminant classifier
for q=1:folds
rand_num = randperm(size(X,1));
X_train = X(rand_num(1:round((1-1/folds)*length(rand_num))),:);
Y_train = Y(rand_num(1:round((1-1/folds)*length(rand_num))),:);

X_test = X(rand_num(round((1-1/folds)*length(rand_num))+1:end),:);
Y_test = Y(rand_num(round((1-1/folds)*length(rand_num))+1:end),:);

LDModels=cell(ms,1);
LSVMModels=cell(ms,1);
QSVMModels=cell(ms,1);
FKNNModels=cell(ms,1);

for j = 1:numel(classes)
    for k=1:round((1-1/folds)*length(rand_num))
        if Y_train(k)==classes(j)
            indx1(k)=1;
        else
            indx1(k)=0;
        end
    end
    indx=indx1';    
    LDModels{j}=fitcdiscr(X_train,indx,'DiscrimType', 'linear','Gamma', 0, 'FillCoeffs', 'off');
    LSVMModels{j}=fitcsvm(X_train,indx,'ClassNames',[false true],'Standardize',true,'KernelFunction','linear');
    QSVMModels{j}=fitcsvm(X_train,indx,'ClassNames',[false true],'Standardize',true,'KernelFunction','polynomial');
    FKNNModels{j}=fitcknn(X_train,indx,'Distance', 'Euclidean','Exponent', [],'NumNeighbors', 5, 'DistanceWeight', 'Equal','Standardize', true);
end

N=size(X_test,1);
Scores=zeros(N,numel(classes));

%linear discriminant classification
for j=1:numel(classes)
    [~,score]=predict(LDModels{j},X_test);
    Scores(:,j)=score(:,2); % Second column contains positive-class scores
end
[~,maxScore]=max(Scores,[],2);

for k=1:size(X_test,1)
    cmp(k)=(maxScore(k)==Y_test(k));
end

accuracy(1,q)=(sum(cmp(:))/length(cmp))*100;

confusionmat=zeros(numberOfClusters);

for i=1:size(X_test,1)
    confusionmat(Y_test(i), maxScore(i))=confusionmat(Y_test(i),maxScore(i))+1;
end

precision(1,q)=(confusionmat(1,1)/(confusionmat(1,1)+confusionmat(1,2)));
recall(1,q)=(confusionmat(1,1)/(confusionmat(1,1)+confusionmat(2,1)));
specificity(1,q)=(confusionmat(2,2)/(confusionmat(2,2)+confusionmat(1,2)));
fscore(1,q)=(2*precision(1,q)*recall(1,q))/(precision(1,q)+recall(1,q));

%Linear SVM classification
for j=1:numel(classes)
    [~,score]=predict(LSVMModels{j},X_test);
    Scores(:,j)=score(:,2); % Second column contains positive-class scores
end
[~,maxScore]=max(Scores,[],2);

for k=1:size(X_test,1)
    cmp(k)=(maxScore(k)==Y_test(k));
end

accuracy(2,q)=(sum(cmp(:))/length(cmp))*100;

confusionmat=zeros(numberOfClusters);

for i=1:size(X_test,1)
    confusionmat(Y_test(i), maxScore(i))=confusionmat(Y_test(i),maxScore(i))+1;
end

precision(2,q)=(confusionmat(1,1)/(confusionmat(1,1)+confusionmat(1,2)));
recall(2,q)=(confusionmat(1,1)/(confusionmat(1,1)+confusionmat(2,1)));
specificity(2,q)=(confusionmat(2,2)/(confusionmat(2,2)+confusionmat(1,2)));
fscore(2,q)=(2*precision(2,q)*recall(2,q))/(precision(2,q)+recall(2,q));


%Quadratic SVM classification
for j=1:numel(classes)
    [~,score]=predict(QSVMModels{j},X_test);
    Scores(:,j)=score(:,2); % Second column contains positive-class scores
end
[~,maxScore]=max(Scores,[],2);

for k=1:size(X_test,1)
    cmp(k)=(maxScore(k)==Y_test(k));
end

accuracy(3,q)=(sum(cmp(:))/length(cmp))*100;

confusionmat=zeros(numberOfClusters);

for i=1:size(X_test,1)
    confusionmat(Y_test(i), maxScore(i))=confusionmat(Y_test(i),maxScore(i))+1;
end

precision(3,q)=(confusionmat(1,1)/(confusionmat(1,1)+confusionmat(1,2)));
recall(3,q)=(confusionmat(1,1)/(confusionmat(1,1)+confusionmat(2,1)));
specificity(3,q)=(confusionmat(2,2)/(confusionmat(2,2)+confusionmat(1,2)));
fscore(3,q)=(2*precision(3,q)*recall(3,q))/(precision(3,q)+recall(3,q));

%Fine KNN classification
for j=1:numel(classes)
    [~,score]=predict(FKNNModels{j},X_test);
    Scores(:,j)=score(:,2); % Second column contains positive-class scores
end
[~,maxScore]=max(Scores,[],2);

for k=1:size(X_test,1)
    cmp(k)=(maxScore(k)==Y_test(k));
end

accuracy(4,q)=(sum(cmp(:))/length(cmp))*100;

confusionmat=zeros(numberOfClusters);

for i=1:size(X_test,1)
    confusionmat(Y_test(i), maxScore(i))=confusionmat(Y_test(i),maxScore(i))+1;
end

precision(4,q)=(confusionmat(1,1)/(confusionmat(1,1)+confusionmat(1,2)));
recall(4,q)=(confusionmat(1,1)/(confusionmat(1,1)+confusionmat(2,1)));
specificity(4,q)=(confusionmat(2,2)/(confusionmat(2,2)+confusionmat(1,2)));
fscore(4,q)=(2*precision(4,q)*recall(4,q))/(precision(4,q)+recall(4,q));
end
classifiers = {"LinearDiscriminant","LinearSVM","QuadraticSVM","FineKNN",'Average'};

Accuracy=[mean(accuracy(1,:)) mean(accuracy(2,:)) mean(accuracy(3,:)) mean(accuracy(4,:))];
avg_acc=mean(Accuracy);
Precision=[mean(precision(1,:)) mean(precision(2,:)) mean(precision(3,:)) mean(precision(4,:))];
avg_precision=mean(Precision);
Recall=[mean(recall(1,:)) mean(recall(2,:)) mean(recall(3,:)) mean(recall(4,:))];
avg_recall=mean(Recall);
Specificity=[mean(specificity(1,:)) mean(specificity(2,:)) mean(specificity(3,:)) mean(specificity(4,:))];
avg_specificity=mean(Specificity);
Fscore=[mean(fscore(1,:)) mean(fscore(2,:)) mean(fscore(3,:)) mean(fscore(4,:))];
avg_fscore= mean(Fscore);

acc = [Accuracy avg_acc];
pre = [Precision avg_precision];
rec = [Recall avg_recall];
spe = [Specificity avg_specificity];
fsco = [Fscore avg_fscore];

s.ClassifierModel=cellstr(classifiers.');
p.ClassifierModel=cellstr(classifiers.');
w.ClassifierModel=cellstr(classifiers.');
e.ClassifierModel=cellstr(classifiers.');
r.ClassifierModel=cellstr(classifiers.');

s.Accuracy=(acc).';
p.Precision=(pre).';
w.Recall=(rec).';
e.Specificity=(spe).';
r.fscore = (fsco).';

ClassificationAccuracyScores = struct2table(s)
ClassificationprecisionScores = struct2table(p)
ClassificationrecallScores = struct2table(w)
ClassificationspecificityScores = struct2table(e)
ClassificationfscoreScores = struct2table(r)

%new sequence classification 
fprintf('New Sequence Classification.... \n');
%uncomment the following lines and replace the ML model(FKNNModels) to LSVMModels/QSVMModels/LDModels 
%in line 276 if required.
%{
[Acnumbers,Z]=newSeqClassify('Test',FKNNModels,feature,numberOfClusters,clusterNames);
m.Genome=cellstr(Acnumbers.');
m.Classified_cluster=cellstr(Z.');
newSeqClassification = struct2table(m)
%}