function [Acnumbers,Z] = newSeqClassify(name1,svm_models,feature,numberOfClusters,clusterNames)
%NEWSEQCLASSIFY Summary of this function goes here
%   Detailed explanation goes here


[AcNumbers1, newSeq] = readFasta1(name1);
totalSeq = length(newSeq);
Acnumbers=AcNumbers1(1:totalSeq);

if feature=="ETC"
    
nmValSH=cell(1,totalSeq);
lg=cell(1,totalSeq);
%Z=cell(1,totalSeq);
newcomplexities=zeros(1,totalSeq);

fprintf('applying ETC for new sequences.... \n');
for j = 1:totalSeq  
    ns = numMappingPP(newSeq{j});
    nmValSH{j} = ns;
    lg{j} = ETC(ns,0);
    newcomplexities(1,j)=lg{j}/(length(newSeq{j})-1);
end
newX=newcomplexities';
newX(:,2)=newX(:,1);
end

if feature=="LZ"

nmValSH=cell(1,totalSeq);%empty cells of 1x148
lg=cell(1,totalSeq);
newcomplexities=zeros(1,totalSeq);

fprintf('applying Lempel-ziv for new sequences .... \n');
for a = 1:totalSeq 
    Seq1{a}=newSeq{a};
    
     [codice, code_bin, code_book]= lempel_ziv(['A' 'T' 'C' 'G'],cell2mat(Seq1(a)));
     lg{a}=length(code_book);
     
     newcomplexities(1,a)=(lg{a}/(length(newSeq{a})))*(log(length(newSeq{a}))/log(4));
end
newX=newcomplexities';
newX(:,2)=newX(:,1);
end

if feature=="ETC+LZ"
    
nmValSH=cell(1,totalSeq);%empty cells of 1x148
lg=cell(1,totalSeq);
newnoretc=zeros(1,totalSeq);

fprintf('applying ETC for new sequences.... \n');
for j = 1:totalSeq  
    ns = numMappingPP(newSeq{j});
    nmValSH{j} = ns;
    lg{j} = ETC(ns,0);
    newnoretc(1,j)=lg{j}/(length(newSeq{j})-1);
end

nmValSH=cell(1,totalSeq);%empty cells of 1x148
lg=cell(1,totalSeq);
newnorlz=zeros(1,totalSeq);

fprintf('applying Lempel-ziv for new sequences .... \n');
for a = 1:totalSeq 
    Seq1{a}=newSeq{a};
    
     [codice, code_bin, code_book]= lempel_ziv(['A' 'T' 'C' 'G'],cell2mat(Seq1(a)));
     lg{a}=length(code_book);
     
     newnorlz(1,a)=(lg{a}/(length(newSeq{a})))*(log(length(newSeq{a}))/log(4));
end

newX=newnoretc';
newX(:,2)=newnorlz';
end


N=size(newX,1);
Scores1=zeros(N,numberOfClusters);
for j=1:numberOfClusters
    [~,score1]=predict(svm_models{j},newX);
    Scores1(:,j)=score1(:,2); % Second column contains positive-class scores
end
[~,maxScore1]=max(Scores1,[],2);

%classification

for k=1:size(newX,1)
    Z(k)=clusterNames(maxScore1(k));
end

end


