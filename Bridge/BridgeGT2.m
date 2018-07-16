data=importdata('Bridge.txt');
 x=data(:,1:2);
 out=data(:,3);
 n=length(x);

disp(n);
p=size(x,2);
Z=x;
m=n;
beta1=100;
beta2=120;

 dis=zeros(n,p,m);
 disbar=zeros(n,m);
 s=zeros(n,m);

for J=1:n
    dis(:,:,J)=abs(x-Z(J,:));
end


mf=SG(x,beta1,beta2);



%Computing dis and disbar
disbar= sum((dis.*mf),2)./sum(mf,2);
               


%Similarity matrix
for i=1:n
    for k=1:m
        if(i~=k)
            s(i,k)=-disbar(i,k);
        end
    end
end


%calculating max and min of similarity matrix

mins=min(min(s));
maxs=s(1,2);
for i=1:n
    for k=1:m
        if(i~=k && s(i,k)>=maxs)
            maxs=s(i,k);
        end
    end
end

            
%preference values
 for i=1:n
   s(i,i)= (mins - 10*(maxs-mins));
end

S=s;


% Initialize messages
N=size(S,1); 
A=zeros(N,N); 
R=zeros(N,N); 
S=(S+(eps*randn(N,N))*(max(S(:))-min(S(:)))); % Remove degeneracies
lam=0.9; % Set damping factor
for i=1:100
    % Compute responsibilities
    Rold=R;
    AS=A+S; 
    [Y,I]=max(AS,[],2);
    for i=1:N 
      AS(i,I(i))=-realmax; 
    end
    [Y2,I2]=max(AS,[],2);
    R=S-repmat(Y,[1,N]);
    for i=1:N 
      R(i,I(i))=S(i,I(i))-Y2(i);
    end
    R=(1-lam)*R+lam*Rold; % Dampen responsibilities

    % Compute availabilities
    Aold=A;
    Rp=max(R,0); 
    for k=1:N 
      Rp(k,k)=R(k,k); 
    end
    A=repmat(sum(Rp,1),[N,1])-Rp;
    dA=diag(A); 
    A=min(A,0); 
    for k=1:N 
      A(k,k)=dA(k);
    end
    A=(1-lam)*A+lam*Aold; % Dampen availabilities
end
E=R+A; % Pseudomarginals
I=find(diag(E)>0); 
K=length(I); % Indices of exemplars
[tmp,c]=max(S(:,I),[],2); 
idx=I(c); % Assignments

uidx=unique(idx);
for it=1:length(uidx)
    idx(uidx(it))= uidx(it);
end


disp([unique(idx),histc(idx,unique(idx))]);
acc=zeros(2,2);
y=zeros(1,n);
ybar=zeros(1,n);

for i=1:n
    if (idx(i)>=1 && idx(i)<=87)
        ybar(i)=3;
    elseif (idx(i)>=88 && idx(i)<=240)
        ybar(i)=6;
    end
end
for i=1:n
    acc(out(i)/3,ybar(i)/3)= acc(out(i)/3,ybar(i)/3)+1;
end

disp(acc);




%%%%%%%%%%% Secondary function %%%%%%%%%%%%%%%%
function mf=SG(x,beta1,beta2)




p=2:0.1:5;
c=(beta1+beta2)/2;
sig=(beta2-beta1)/7;
syms m;
gauss(m)=gaussmf(m,[sig c]);
boss=matlabFunction(gauss);
gaussinv=finverse(gauss);
F=matlabFunction(gaussinv);
L=zeros(5,2);
for I=1:5
        L(I,1)=beta1+beta2-F(I/5);
        L(I,2)=F(I/5);
end


n=length(x);
p=size(x,2);
Z=x;
m=n;


beta=L;
u=zeros(5,1);
for I=1:5
    u(I,1)=boss(F(I/5));
end

points=zeros(5,n,p,n);

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
for l=1:1
    
sum21=sum((dis1.*prevmf1),2);
sum31=sum(prevmf1,2);
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

     points(iter,:,:,:)= ((mf2+mf1)/2)*u(iter);


 %%%%%%%%%%%%%%%%%%%major loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
end

mf=sum(points,1);
mf=reshape(mf,n,p,n);


end