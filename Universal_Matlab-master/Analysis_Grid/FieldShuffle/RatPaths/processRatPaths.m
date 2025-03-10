%Script to process the raw rat paths in this folder. In each case reads in
%the list of .mat files which should contain full tint data sets for the
%pos strucutre (i.e. pos.xy). Then normalises this so both x and y lie
%between 0 and 100 (i.e. ppm is set to 100). NB. all this assumes the paths
%are 1m. Then cuts the paths into 1min chunks and concatanates into a
%single large 3D mat (single precision to save space).

chunkLgt        =50*60; %Chunk length in pos samp at 50hz


fileNm          =dir;
fileNm          =fileNm(3:12);

posChunks       =single([]);
for nn          =1:length(fileNm)
   
    load(fileNm(nn).name); %loads pos.xy
    
    %Convert to be between 0 and 100
    pos.xy          =bsxfun(@minus, pos.xy, min(pos.xy));
    pos.xy          =bsxfun(@rdivide, pos.xy, max(pos.xy)) .*100;
 
    %Now chunk into 1min (50*60) pos sample chunks
    xy              =pos.xy;
    clear pos
    nChunk          =floor(length(xy)/chunkLgt);
    xy              =xy(1:chunkLgt*nChunk,:);
    xy              =reshape(xy, [chunkLgt,2, nChunk]);
    posChunks       =cat(3, posChunks,single(xy));
    clear xy nChunk
end

%Finally save posChunks which is a chunkLgt x 2 x nChunk 3D mat containing
%normalised xy pairs of pos data for 1minute (nomrally) of behaviour
save posChunks posChunks