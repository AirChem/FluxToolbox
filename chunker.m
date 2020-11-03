% function chunks = chunker(index)
%input is a vector containing indices for all data chunks of a certain type (e.g. cals or zeroes)
%output is a 2-column matrix containing start and stop indices for each chunk
%Edited from the MBO2006 function autoindexer.m.
%070707 GMW

function chunks = chunker(index)

if isempty(index)
    chunks = [];
else
    index = index(:);
    chunks = [];%this matrix will contain the indices
    j=find(diff(index)~=1);
    chunkstart = [index(1);index(j+1)];
    chunkstop = [index(j);index(end)];
    chunks = [chunkstart chunkstop];
end