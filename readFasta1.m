function [AcNmb, Seq] = readFasta1(dataSet)
% read fasta files in folder DataBase/'dataSet'
       
    path = pwd;
    dbPath = strcat(path,'\','DataBase');
    cd(dbPath);
   
        dbPath = strcat(dbPath,'\',dataSet);
        folderInfo = dir (dbPath);
        Seq = {};
        AcNmb = {};    

        if(isempty(folderInfo))
             error('Error : DataSet does not exist.')
        else
            cd(dbPath); 

            folderInfo = dir (dbPath);
            %for i=3:length(folderInfo)
                cd(dbPath);                
                pts = length(folderInfo)-2;                
                seqTemp = cell(1,pts);
                acTemp = cell(1,pts);
                for j=3:length(folderInfo)
                    [Header, Sequence] = fastaread(folderInfo(j).name);    
                    Sequence = regexprep(Sequence,'[^A,^C, ^G, ^T]','','ignorecase');
                    seqTemp{j-2} = Sequence;
                    acTemp{j-2} = Header; 
                end
                Seq = [Seq seqTemp];
                AcNmb = [AcNmb acTemp]; 
       end
        cd(path);
end

