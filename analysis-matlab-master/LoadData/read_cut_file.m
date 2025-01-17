function [exact, unique_clusters, n_user_clusters] = read_cut_file(flnm)
    
if ~isempty(regexp(flnm,'.clu.[0-9]+$', 'once'))
    [exact, n_user_clusters] = read_cut_file_clu(flnm);
elseif ~isempty(regexp(flnm,'.cut$', 'once'))
    [exact, n_user_clusters] = read_cut_file_cut(flnm);
end

unique_clusters = unique(exact);
end

function [exact, n_user_clusters] = read_cut_file_clu(flnm)
    fid = fopen(flnm);
    n_user_clusters = textscan(fid,'%f',1);
    n_user_clusters = n_user_clusters{:};
    exact = textscan(fid,'%f',Inf);
    exact = exact{:}-1; % -1 to match what tint shows
end

function [exact, n_user_clusters] = read_cut_file_cut(flnm)
% Only the exact cut is extracted at present

TOP_FORMAT = 'n_clusters: %f n_channels: %f n_params: %f times_used_in_Vt: %f %f %f %f';
CLUSTER_FORMAT = ' cluster: %f center: %f %f %f %f %f %f %f %f  min: %f %f %f %f %f %f %f %f max: %f %f %f %f %f %f %f %f';
START_REGEXP = 'spikes: ([0-9]*)'; %regexp here rather than textscan to allow us to jump over arbitrary trial name

fid = fopen(flnm);

if fid ~= -1
    %get the stuff at the top, mostly ignored
    top = textscan(fid,TOP_FORMAT,'Whitespace','\t\n\r ','MultipleDelimsAsOne',1);
    [n_clusters, n_channels,n_params,vt1 vt2 vt3 vt4] = top{:};
    clusters = textscan(fid,CLUSTER_FORMAT,'Whitespace','\t\n\r ','MultipleDelimsAsOne',n_clusters);
    n_spikes = regexp(fgets(fid),START_REGEXP,'tokens'); 
    n_spikes = str2double(n_spikes{1});
    
    %get the cut
    exact = textscan(fid,'%f',n_spikes);
    exact = exact{:};
    
    fclose(fid);
    n_user_clusters = n_clusters;
end

end