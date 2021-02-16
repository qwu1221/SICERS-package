function [signode,GLPstat,P_value,Tmax_5prct,max_T0]=PermTest_pval(Wt,Wt2,Cindx,CID,thresh,M,r)
%% This function is to test the significance of covariate-related subnetworks


%%% Inputs:

%%%%%  Wt,Wt2:  Network data for patients and health controls
%%%%%  Cindx:   The cluster/subnetwork index of every non-isolated node
%%%%%  CID:     The cluster/subnetwork index of every cluster in a power descending
%%%%%           order. i.e. CID(1) will be the cluster index of the most
%%%%%           concentrated cluster
%%%%%  thresh:  A threshold on the p-values
%%%%%  M:       The number of label permutation
%%%%%  r:       The cutoff to binarize the networks


%%% Outputs:

%%%%%  signode: Nodes in significant subnetworks
%%%%%  GLPstat: Test statistic for each subnetwork
%%%%%  P_value: The p value for each subnetwork
%%%%%  Tmax_5prct: The 95 percentage for test statistic under null
%%%%%  max_T0:  The largest test statistic of subnetworks in the network
%%%%%  data

n1=size(Wt,2);
n2=size(Wt2,2);

[h1,p1]=ttest2(Wt',Wt2');
Wt_c=[Wt Wt2];
W1=squareform(-log(p1));
Cindx0=Cindx;
CID0=CID;
K_0=size(CID,2); %number of clusters
allinx=1:size(W1,1);

phat=sum(-log(p1)>r)/size(p1,2);

T_k0=zeros(K_0,1);
qhat_k0=zeros(K_0,1);
k_0 = zeros(K_0,1);% number of nodes in kth cluster
for c = 1:K_0
    idx_c = find(Cindx==c);
    nidx_c = size(idx_c,2);
    k_0(c) = nidx_c;
    nedges = nidx_c*(nidx_c-1)/2;
    G_k = W1(idx_c,idx_c);
    if(nidx_c<=4)
        T_k0(c) = 0;
    else
        qhat_k0(c) = sum(squareform(G_k)>r)/size(squareform(G_k),2);
        T_k0(c)=nidx_c^2/(2/(qhat_k0(c)-phat)^2+2/3/(qhat_k0(c)-phat));
    end
    close all;
end


T_max = zeros(M,1);
error_idx = zeros(M,1);
for m=1:M
    try
    ar=randperm(n1+n2);
    Wt_c=Wt_c(:,ar);
    Wt1_c=Wt_c(:,1:n1);
    Wt2_c=Wt_c(:,(1+n1):(n1+n2));
    [h1,p1]=ttest(Wt1_c',Wt2_c');
    
    W1=squareform(-log(p1));
    phat = sum(-log(p1)>r)/size(p1,2);
    
    nlogp=squareform(W1);
    [Cindx,CID,Clist]=SICERS(nlogp,thresh,0,3);
    K=size(CID,2); %number of clusters
    G=squareform(nlogp);
    T_k=zeros(K,1);
    qhat_k=zeros(K,1);;
    for c=1:K
        idx_c=find(Cindx==c);
        nidx_c=size(idx_c,2);
        nedges=nidx_c*(nidx_c-1)/2;
        G_k=G(idx_c,idx_c);
        if(nidx_c<=4)
            T_k(c)=0;
        else
             qhat_k(c) = sum(squareform(G_k)>r)/size(squareform(G_k),2);
             T_k(c)=nidx_c^2/(2/(qhat_k(c)-phat)^2+2/3/(qhat_k(c)-phat));
        end
        close all;
    end
    T_max(m)=max(T_k);
    catch 
        error_idx(m)=1;
    end
end
T_max=T_max(error_idx==0);
Tmax_5prct = prctile(T_max,95);
max_T0 = max(T_k0);
[diff_value,diff_ind]=min(abs(sort(T_max)-max(T_k0)));
P_value=1-diff_ind/M;


%% account for multiple significant subnetworks
signode={};
SCinx = find(T_k0>Tmax_5prct);
noc = size(SCinx,1);
GLPstat=zeros(noc,1);
for j=1:noc
    signode{j}=allinx(ismember(Cindx0,SCinx(j)));
    GLPstat(j) = T_k0(SCinx(j));
end
end