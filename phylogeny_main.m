%Phylogenetic tree creation using ETC/LZ distance matrix.
%clear workspace
close all;
clear all;
clc ;
set(0, 'DefaultFigureRenderer', 'painters')

%read fasta files
dataSet = 'Mammals';
selectedFolder = dataSet;
fprintf('Reading sequences .... \n');
[AcNmb, Seq, numberOfClusters, clusterNames, pointsPerCluster] = readFasta(dataSet);
totalSeq = length(Seq);

%calculate length stats
[maxLen, minLen, meanLen, medLen] = lengthCalc(Seq);

nmValSH=cell(1,totalSeq);
lg=cell(1,totalSeq);
noretc=cell(1,totalSeq);

feature=input("Enter complexity measure(ETC/LZ): ");

if (feature=="ETC")
fprintf('Generating numerical sequences, applying ETC .... \n');
parfor a = 1:totalSeq  %parallel for loop
    ns = numMappingPP(Seq{a});
    nmValSH{a} = ns;
    lg{a} = ETC(ns,0);
    noretc{a}=lg{a}/(length(Seq{a})-1);
end

%distance calculation
fprintf('Computing Distance matrix .... \n');
disMat = zeros(totalSeq);
comETC = zeros(totalSeq);
for i=1:totalSeq
    for j=i:totalSeq
        c1=[cell2mat(nmValSH(i)), cell2mat(nmValSH(j))];
        c2=[cell2mat(nmValSH(j)), cell2mat(nmValSH(i))];
        c1ETC=ETC(c1,0)/(length(c1)-1);
        c2ETC=ETC(c2,0)/(length(c2)-1);
        comETC(i,j)=c1ETC;
        comETC(j,i)=c2ETC;
        disMat(i,j)= (c1ETC-cell2mat(noretc(i))+c2ETC-cell2mat(noretc(j)))/ (0.5*(c1ETC +c2ETC));
        disMat(j,i)=disMat(i,j);
    end
end

elseif (feature=="LZ")
lg=cell(1,totalSeq);

fprintf('applying Lempel-ziv for sequences .... \n');
for a = 1:totalSeq  %parallel for loop
    
    [codice, code_bin, code_book]= lempel_ziv(['A' 'T' 'C' 'G'],cell2mat(Seq(a)));
    lg{a}=length(code_book);
     
    norlem{a}=(lg{a}/(length(Seq{a})))*(log(length(Seq{a}))/log(4));
       
end

%distance calculation
fprintf('Computing Distance matrix .... \n');
disMat = zeros(totalSeq);
comlem = zeros(totalSeq);
for i=1:totalSeq
    for j=i:totalSeq
        c1=[cell2mat(Seq(i)), cell2mat(Seq(j))];
        c2=[cell2mat(Seq(j)), cell2mat(Seq(i))];
        [codice, code_bin, code_book]=lempel_ziv(['A' 'T' 'C' 'G'],c1);
        c1ziv=(length(code_book)/(length(c1)))*(log(length(c1))/log(4));
        [codice, code_bin, code_book]=lempel_ziv(['A' 'T' 'C' 'G'],c2);
        c2ziv=(length(code_book)/(length(c2)))*(log(length(c2))/log(4));
        comlem(i,j)=c1ziv;
        comlem(j,i)=c2ziv;
        disMat(i,j)= (c1ziv-cell2mat(norlem(i))+c2ziv-cell2mat(norlem(j)))/ (0.5*(c1ziv +c2ziv)) ;
        disMat(j,i)= disMat(i,j);
    end
end

else
    fprintf('Wrong feature provided! \n');
end

minval=min(disMat,[],'all');
for i=1:totalSeq
    for j=1:totalSeq
        disMat(i,j)=(disMat(i,j)-minval);
    end
end

maxval=max(disMat,[],'all');
for i=1:totalSeq
    for j=1:totalSeq
        disMat(i,j)=(disMat(i,j)/maxval);
    end
end

for i=1:totalSeq
    disMat(i,i)=0;
end




%Phylogenetic Tree
% plotted for datasets with atmost 60 sequences
if(totalSeq<=60)
    fprintf('Creating Phylogenetic Tree .... \n');
    UPGMAtree = seqlinkage(disMat,'UPGMA',AcNmb);
    plot(UPGMAtree);
end






