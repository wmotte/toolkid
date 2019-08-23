% Script that generates indices and tables for DKI

%% First create the index table
index_table=zeros(81,4);

lin_idx=0;
for i_idx=1:3
    for j_idx=1:3
        for k_idx=1:3
            for l_idx=1:3
                lin_idx=lin_idx+1;
                index_table(lin_idx,:)=[i_idx j_idx k_idx l_idx];
            end
        end
    end
end

%% Next, from the fact that we know the order in which DKI parameters are
%% calculated, calculate reverse indices into the table

% kvec is of the form k1111 k2222 k3333 k1112 k1113 k1222 k2223 k1333 k2333
% k1122 k1133 k2233 k1123 k1223 k1233
kvec_indices=[1 1 1 1; 2 2 2 2; 3 3 3 3;1 1 1 2; 1 1 1 3; 1 2 2 2; 2 2 2 3; 1 3 3 3; 2 3 3 3; 1 1 2 2; 1 1 3 3; 2 2 3 3; 1 1 2 3; 1 2 2 3; 1 2 3 3];


sorted_index_table=sort(index_table,2); % reverse rows
sorted_kvec_indices=sort(kvec_indices,2); % reverse rows

index_table_str=int2str(sorted_index_table);
kvec_indices_str=int2str(sorted_kvec_indices);

% Get mapping from kvec_indices to index_table
kvec_to_table=zeros(81,1);
for idx=1:size(kvec_indices_str,1)
    matches=strcmp(kvec_indices_str(idx,:),cellstr(index_table_str));
%     nnz(matches)
    kvec_to_table(matches)=idx;
end


% 

indices_wt1111=[sub2ind([3 3],index_table(:,1),repmat(1,[81 1])) sub2ind([3 3],index_table(:,2),repmat(1,[81 1])) sub2ind([3 3],index_table(:,3),repmat(1,[81 1]))  sub2ind([3 3],index_table(:,4),repmat(1,[81 1]))];
indices_wt2222=[sub2ind([3 3],index_table(:,1),repmat(2,[81 1])) sub2ind([3 3],index_table(:,2),repmat(2,[81 1])) sub2ind([3 3],index_table(:,3),repmat(2,[81 1]))  sub2ind([3 3],index_table(:,4),repmat(2,[81 1]))];
indices_wt3333=[sub2ind([3 3],index_table(:,1),repmat(3,[81 1])) sub2ind([3 3],index_table(:,2),repmat(3,[81 1])) sub2ind([3 3],index_table(:,3),repmat(3,[81 1]))  sub2ind([3 3],index_table(:,4),repmat(3,[81 1]))];
indices_wt2233=[sub2ind([3 3],index_table(:,1),repmat(2,[81 1])) sub2ind([3 3],index_table(:,2),repmat(2,[81 1])) sub2ind([3 3],index_table(:,3),repmat(3,[81 1]))  sub2ind([3 3],index_table(:,4),repmat(3,[81 1]))];
indices_wt1133=[sub2ind([3 3],index_table(:,1),repmat(1,[81 1])) sub2ind([3 3],index_table(:,2),repmat(1,[81 1])) sub2ind([3 3],index_table(:,3),repmat(3,[81 1]))  sub2ind([3 3],index_table(:,4),repmat(3,[81 1]))];
indices_wt1122=[sub2ind([3 3],index_table(:,1),repmat(1,[81 1])) sub2ind([3 3],index_table(:,2),repmat(1,[81 1])) sub2ind([3 3],index_table(:,3),repmat(2,[81 1]))  sub2ind([3 3],index_table(:,4),repmat(2,[81 1]))];

dki_tabs.index_table=index_table;
dki_tabs.kvec_to_table=kvec_to_table;
dki_tabs.indices_wt1111=indices_wt1111;
dki_tabs.indices_wt2222=indices_wt2222;
dki_tabs.indices_wt3333=indices_wt3333;
dki_tabs.indices_wt2233=indices_wt2233;
dki_tabs.indices_wt1133=indices_wt1133;
dki_tabs.indices_wt1122=indices_wt1122;


save('dki_table_indices.mat','dki_tabs');

