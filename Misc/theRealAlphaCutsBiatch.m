

data=split(importdata('iris dataset.txt'),",");
 size(data);
 X=data(:,1:4);
 Y=data(:,5);
 x=str2double(X);
beta1=2;
beta2=5;
p=2:0.1:5;
c=(beta1+beta2)/2;
sig=(beta2-beta1)/7;
syms m;
gauss(m)=gaussmf(m,[sig c]);%this forms a guass function in the form of a string
boss=matlabFunction(gauss);%forms the actual usable gauss function
gaussinv=finverse(gauss);
F=matlabFunction(gaussinv);%inverse gauss function
L=zeros(5,2);
update=0.1;
for I=1:5
        L(I,1)=beta1+beta2-F(I/5);
        L(I,2)=F(I/5);
end


n=length(x);
p=size(x,2);
Z=x;
m=n;

mfup=rand(n,p,n);
mfdown=rand(n,p,n);

beta=L;
u=zeros(5,1);
for I=1:5
    u(I,1)=boss(F(I/5));
end

points=zeros(10,n,p,n);
alpha=zeros(10,n,p,n);

for iter=1:5
%%%%%%%%%%%%%%%%%%%%%major loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Initalizations for membership function corresponding to beta1
mf1=ones(n,p,m);
prevmf1=ones(n,p,m);
dis1=zeros(n,p,m);
disbar1=zeros(n,p,m);
devbar1=zeros(n,m);
dev1=zeros(n,p,m);
prevdisbar1=zeros(n,m);
comparison2=1;

%Initializations for membership function corresponding to beta2
mf2=ones(n,p,m);
prevmf2=ones(n,p,m);

disbar2=zeros(n,p,m);
devbar2=zeros(n,m);
dev2=zeros(n,p,m);

comparison4=1;

%initializing dis and disbar

for J=1:n
    dis1(:,:,J)=abs(x-Z(J,:));
end
        disbar1=squeeze(sum(dis1,2)/p);
        prevdisbar2=prevdisbar1;

%iterations
for l=1:20
    
sum21=sum((dis1.*prevmf1),2);
sum31=sum(prevmf1,2);
sum31=(sum31<=0)+sum31;
disbar1=sum21./sum31;
disbar1(isnan(disbar1))=0;

sum22=sum((dis1.*prevmf2),2);
sum32=sum(prevmf2,2);
disbar2=sum22./sum32;
disbar2(isnan(disbar2))=0;

B1=repmat(disbar1,1,p,1);
B2=repmat(disbar2,1,p,1);

dev1=abs(dis1-B1);
dev2=abs(dis1-B2);

devbar1=(sum(dev1,2))/p;
devbar2=(sum(dev2,2))/p;

P1=repmat(devbar1,1,p,1);
P2=repmat(devbar2,1,p,1);

    mf1=exp(1).^(-(dev1.^beta(iter,1))./P1.^beta(iter,1));
    mf2=exp(1).^(-(dev2.^beta(iter,2))./P2.^beta(iter,2));
    mf1(isnan(mf1))=0;
    mf2(isnan(mf2))=0;


prevmf1=mf1;
prevmf2=mf2;
comparison2=0;

G1=abs(prevdisbar1-disbar1);
G2=abs(prevdisbar2-disbar2);

if sum(sum(sum(G1>0.01)))~=0 || sum(sum(sum(G2>0.01)))~=0
    comparison2=1;
end

prevdisbar1=disbar1;
prevdisbar2=disbar2;

end

     points(iter,:,:,:)= mf1;
     points(11-iter,:,:,:)= mf2;
     
    alpha(iter,:,:,:)=repmat(u(iter),n,p,n);
    alpha(11-iter,:,:,:)=repmat(u(iter),n,p,n);

 %%%%%%%%%%%%%%%%%%%major loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
end

points(5,:,:,:)=[];
alpha(5,:,:,:)=[];
points=(~(points<0.001)).*points;

secgrade=zeros(n,p,n); 
mftemp=mfdown;
sum5=zeros(n,p,n);
sum6=sum5;
parfor i=1:n
    tic
    for k=1:n
        if i~=k
        for j=1:p
%           if max(points(:,i,j,k))==min(points(:,i,j,k))
%               sum5(i,j,k)=0;b 
%               sum6(i,j,k)=0;
%           else
          if round(max(points(:,i,j,k))*10000)/10000~=round(min(points(:,i,j,k))*10000)/10000
              secgrade= fit(points(:,i,j,k),alpha(:,i,j,k),'gauss3');
       
              while mfdown(i,j,k)<=mftemp(i,j,k) && mftemp(i,j,k)<=mfup(i,j,k)
               
                 sum5(i,j,k)= sum5(i,j,k)+ mftemp(i,j,k)*abs(secgrade(mftemp(i,j,k)));
                 sum6(i,j,k)= sum6(i,j,k) + abs(secgrade(mftemp(i,j,k)));
                 mftemp(i,j,k)=mftemp(i,j,k)+update;
                
              end
          else
              sum5(i,j,k)=1;
              sum6(i,j,k)=1;
          end
        end
        end
    end
    toc
end

% 
% mftemp=mfdown;
% 
% sum5=0;
% sum6=0;
% 
% if count~=9
%      while mfdown<=mftemp && mftemp<=mfup
%         sum5= sum5+ mftemp*secgrade(mftemp);
%         sum6= sum6 + secgrade(mftemp);
%         mftemp=mftemp+update;
%       end
% else
%     sum5=1;
%     sum6=1;
% end
% if sum5==0 && sum6==0
%         sum6=1;
%     end
% toc
time=toc;