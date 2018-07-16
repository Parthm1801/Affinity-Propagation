
X=getdata('s1_entro.jpg','s1_med.jpg','s1_gauss.jpg','s1_mask.jpg',5,3);
x=X(:,1:3);
out=X(:,4);
n=length(x);
p=size(x,2);
Z=x;
m=n;
mf=ones(n,p,n);
prevmf1=ones(n,p,n);
dis=zeros(n,p,m);
disbar=zeros(n,m);
devbar=zeros(n,m);
dev=zeros(n,p,m);
beta1=0.1;
beta2=1;
prevdisbar=zeros(n,m);
comparison2=1;
s=zeros(n,m);

%initializing prevmf
for i=1:n
        prevmf1(i,:,i)= 0;
  
end

mf=prevmf1;
prevmf2=prevmf1;
%initializing dis and disbar
for J=1:n
    dis(:,:,J)=abs(x-Z(J,:));
end
        disbar=squeeze(sum(dis,2)/p);


%iterations
for it=1:4
    
sum2=sum((dis.*mf),2);
sum3=sum(mf,2);
disbar=sum2./sum3;


%Dev and devbar
B=repmat(disbar,1,p,1);
dev=abs(dis-B);
devbar=(sum(dev,2))/p;

P=repmat(devbar,1,p,1);

%membershipfunction1
prevmf1=exp(1).^(-(dev.^beta1)./(P.^beta1));
prevmf1(isnan(prevmf1))=0;
             
           

%membershipfunction2
prevmf2=exp(1).^(-(dev.^beta2)./(P.^beta2));
prevmf2(isnan(prevmf2))=0;


mf=(prevmf1+prevmf2)/2;
% comparison2=0;
% 
% G=abs(prevdisbar-disbar);

% if sum(G(:)>0.01)~=0
%     comparison2 =1;
% end

prevdisbar=disbar;


end

    

%Similarity matrix
s=-squeeze(disbar);

for i=1:n
    s(i,i)=0;
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
%
% I=(mean(mean(s))*eye(n))/4;
% s=s+I;
 for i=1:n
   s(i,i)= (mins - 60*(maxs-mins));
end

S=s;

% Initialize messages
N=size(S,1); 
A=zeros(N,N); 
R=zeros(N,N); 
S=(S+(eps*randn(N,N))*(max(S(:))-min(S(:)))); % Remove degeneracies
lam=0.9; % Set damping factor
for it=1:100
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
% disp(idx)
uidx=unique(idx);

for i=1:length(uidx)
    idx(uidx(i))=uidx(i);
end

disp([unique(idx),histc(idx,unique(idx))]);

%mkplot(x,idx);

acc=zeros(3);
y=zeros(1,n);
ybar=zeros(1,n);

for i=1:n
    if (idx(i)>=1 && idx(i)<=130)
        ybar(i)=0;
    elseif (idx(i)>=131 && idx(i)<=1032)
        ybar(i)=127.5;
    elseif (idx(i)>=1033 && idx(i)<=1521)
        ybar(i)=255;

    end
end
for i=1:n
    acc((out(i)/127.5)+1,(ybar(i)/127.5)+1)= acc((out(i)/127.5)+1,(ybar(i)/127.5)+1)+1;
end

disp(acc);
mkplot(X,idx)

function x=getdata(filename2,filename3,filename4,filename1,interval,n)

% I=imread(filename1);
% I=double(I);
% miny=min(min(I))*ones(size(I));
% maxy=max(max(I))*ones(size(I));
% img=(I-miny)./(maxy-miny);
% map=zeros(size(I));
% for i=1:n
% 
%     map=map+((i-1)*255/(n-1))*(img(:,:,1)>=((i-1)/(n-1))-(1/(2*(n-1))) & img(:,:,1)<((i-1)/(n-1))+(1/(2*(n-1))));
% end
map=imread(filename1);
  map=double(map);
 map=(map>100 & map<180)*127.5 + (map<100)*0 + (map>180)*255;
 % map=(map<100)*0 + (map>100)*255;
check=size(map,1);

indarr=[];
val=zeros(n,1);
 for row=1:interval:floor(size(map,1)-interval)
        for col=1:interval:floor(size(map,2)-interval)
            i=row+floor(rand*interval);
            j=col+floor(rand*interval);
            temp=map(i,j);
            
            for it=1:n
               if temp==((it-1)*255/(n-1)) 
                temparr=(j-1)*check+i;
                indarr(end+1)=temparr;
                val(it,1)=val(it,1)+1;
               end 
            end

        end
 end
 
 
 I1=getfeatures(filename2,indarr);
 I2=getfeatures(filename3,indarr);
 I3=getfeatures(filename4,indarr);
 Out=map(indarr(1,:));
  Out=Out';
 
 x=[I1 I2 I3 Out];
 [dummy index]=sort(x(:,4));
 x=x(index,:);
 
end


 function I=getfeatures(filename,arr)
  i1=imread(filename);
  i1=double(i1);
  I=i1(arr(1,:));
  I=I';
 end
 
 function mkplot(X,idx)
uidx=unique(idx);
sidx=zeros(size(idx,1),size(idx,2));
% for i=1:length(uidx)
%     sidx=sidx+i*(idx==uidx(i));
% end
color=idx/sum(uidx);
% color=sidx/max(max(sidx));
scatter3(nrmlz(X(:,1)),nrmlz(X(:,2)),nrmlz(X(:,3)),100,color,'filled');
 end
 
 
 function Nx=nrmlz(x)
    Nx=(x-min(x))./(max(x)-min(x));
 end

 