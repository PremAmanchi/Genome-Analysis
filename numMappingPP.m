function [ numSeq ] = numMappingPP( sq )
%function for Purine/Pyramidine representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Kanike Jayanth  %
% Department of Electronics and Communication,%
% Amrita Vishwa Vidyapeetham %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%purines A,G=-1 ; pyramidines C,T=1
    len = length(sq);  
    numSeq = zeros(1,len,'double');
    for K = 1:len
       t = sq(K);
       if(strcmpi(t,'A'))
           numSeq(K) = -1;
       elseif(strcmpi(t,'C'))
           numSeq(K) = 1;
       elseif(strcmpi(t,'G'))
           numSeq(K) = -1; 
       elseif(strcmpi(t,'T'))
           numSeq(K) = 1;
       end           
   end   
end

