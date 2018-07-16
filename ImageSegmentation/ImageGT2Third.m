[x indices]=get_image_data_v3('s3.jpg');
 %x=importdata('x2.txt');
 %indices=importdata('indices3.txt');
 n=length(x);
p=size(x,2);
Z=x;
m=n;
beta1=0.1;
beta2=1;

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
   s(i,i)= (mins - 6*(maxs-mins));
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
for i=1:length(uidx)
    idx(uidx(i))=uidx(i);
end


disp([unique(idx),histc(idx,unique(idx))]);

show_image_v3(x,idx,uidx,'s3.jpg',indices)
%mkplot(x,idx);


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
sumu=0;
points=zeros(5,n,p,n);

for iter=1:5
%%%%%%%%%%%%%%%%%%%%%major loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Initalizations for membership function corresponding to beta1
mf1=ones(n,p,m);
prevmf1=ones(n,p,m);
dis3=zeros(n,p,m);
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
    dis3(:,:,J)=abs(x-Z(J,:));
end
        disbar1=squeeze(sum(dis3,2)/p);
        prevdisbar2=prevdisbar1;

%iterations
for l=1:20
    
sum21=sum((dis3.*prevmf1),2);
sum31=sum(prevmf1,2);
disbar1=sum21./sum31;
disbar1(isnan(disbar1))=0;

sum22=sum((dis3.*prevmf2),2);
sum32=sum(prevmf2,2);
disbar2=sum22./sum32;
disbar2(isnan(disbar2))=0;

B1=repmat(disbar1,1,p,1);
B2=repmat(disbar2,1,p,1);

dev1=abs(dis3-B1);
dev2=abs(dis3-B2);

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

     points(iter,:,:,:)= ((mf2+mf1)/2)*u(iter,1);
    sumu=sumu+u(iter,1);

 %%%%%%%%%%%%%%%%%%%major loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
end

mf=sum(points,1)/sumu;
mf=reshape(mf,n,p,n);


end


function [X all_indices]=get_image_data_v3(filename);
   
    I=imread(filename);
    if size(I,3)==3
        I=rgb2gray(I);
    end
    I=double(I);
    I=(I-min(min(I)))./(max(max(I))-min(min(I)));
    X=[];
    I1=double(imread('s3_entro.jpg'));%gaussian blur
    I2=double(imread('s3_med.jpg'));%median filtering
    I3=double(imread('s3_gauss.jpg'));%entropy
    I_mask=double(imread('s3_mask.jpg'));
    I_mask(find(125>I_mask))=0;
    I_mask(find(I_mask>200))=1;
    I_mask(find(((200>=I_mask).*I_mask)>=125))=0.5;
    all_indices=[];
    for unique_pixel=unique(I_mask)'
        X_temp=[];
        indices=find(I_mask==unique_pixel);
        indices_2=ceil(rand(1,100)*size(indices,1));
        indices_2=indices(indices_2);
        all_indices=[all_indices; indices_2];
        X_temp=[X_temp I1(indices_2)]    ;
        X_temp=[X_temp I2(indices_2)]    ;
        X_temp=[X_temp I3(indices_2)]    ;
        X=[X ;X_temp];
    end 
    %X=X';
    
   
   %scatter3(X_norm(:,1),X_norm(:,2),X_norm(:,3));
end

function show_image_v3(X,idx,uidx,filename,indices)
    color=zeros(size(idx));
    for i=1:length(uidx)
        color=color+(idx==uidx(i))*((i-1)/(length(uidx)-1));
    end
    color=color*255;
    I_mask=double(imread('s3_mask.jpg'));
    I_mask(find(125>I_mask))=0;
    I_mask(find(I_mask>200))=1;
    I_mask(find(((200>=I_mask).*I_mask)>=125))=0.5;
    I_mask=uint8(I_mask*255);
    max_classifications=0;
    unique_color_all=unique(color);
    negative_indices=0;
    while size(unique_color_all)<3
        negative_indices=negative_indices-1;
        unique_color_all=[unique_color_all; negative_indices];
    end
    for unique_color=perms(unique_color_all)'
        %correct_classifications=sum(color(1:100)==unique_color(1))+sum(color(101:200)==unique_color(2))+sum(color(201:300)==unique_color(3));
        correct_classifications=sum(color(1:100)==unique_color(1))+sum(color(101:200)==unique_color(2));
        if max_classifications<correct_classifications
            max_classifications=correct_classifications;
        end
    end
    error_rate=(1-max_classifications/200)*100
    X=X*255;
create_image_2(X,color,indices,I_mask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_whole=uint8(create_image_2(X,color,indices,I_mask));

    indices_list=cell(size(unique_color_all));
    unique_color=unique(X_whole)';
    for unique_color_index=1:size(unique_color,2)
        indices=find(X_whole==unique_color(unique_color_index));
        indices_list(unique_color_index)={indices};
    end
    max_classifications=0;
    confusion_matrix_2=zeros(size(unique_color,1),size(unique_color,1));
    if sum(unique_color_all>=0)
        unique_color_all=zeros(2,1)
        unique_color_all(1)=0;
        unique_color_all(2)=128;
        %unique_color_all(3)=255;
    end
    for unique_color=perms(uint8(unique_color_all))'
        correct_classifications=0;
        for unique_color_index=1:size(unique_color,1)
            X_whole(cell2mat(indices_list(unique_color_index)))=unique_color(unique_color_index);
        end
        correct_classifications=sum(sum(X_whole==I_mask));
        if max_classifications<correct_classifications
            max_classifications=correct_classifications;
            X_temp_whole=X_whole;
        end
    end
    X_whole=X_temp_whole;
    error_rate_whole=(1-max_classifications/(size(I_mask,1)*size(I_mask,2)))*100
   confusion_matrix=zeros(2); 
       unique_color=unique(X_whole);
      for row_index=1:size(confusion_matrix,1)
          for column_index=1:size(confusion_matrix,2)
              confusion_matrix_2(row_index,column_index)=sum(sum(((X_whole==unique_color(column_index))*2-1)==(I_mask==unique_color(row_index))));             
          end
      end
      confusion_matrix_2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function mkplot(X,idx)
uidx=unique(idx);
sidx=zeros(size(idx,1),size(idx,2));
% for i=1:length(uidx)
%     sidx=sidx+i*(idx==uidx(i));
% end
color=idx/sum(uidx);
% color=sidx/max(max(sidx));

hold on
plotcube([1 1 1],[0 0 0],0);
scatter3(nrmlz(X(:,1)),nrmlz(X(:,2)),nrmlz(X(:,3)),100,color,'f');
xlabel('Entropy','FontSize',14,'FontWeight','bold');
ylabel('Median Filtered Intensity','FontSize',14,'FontWeight','bold');
zlabel('Gaussian Blur','FontSize',14,'FontWeight','bold');hold off
end
 
 
 function Nx=nrmlz(x)
    Nx=(x-min(x))./(max(x)-min(x));
 end
 
 

function X_whole=create_image_2(X,color,indices,I_mask)
X_whole=zeros(size(I_mask));    
    for row_index=1:size(I_mask,1)
        for column_index=1:size(I_mask,2)
            min_distance=inf;
            min_point=zeros(9,1);
            flag=zeros(size(indices,1),1);
            for iter=1:9
                for index=1:size(indices,1)
                    column=floor((indices(index)-0.1)/size(I_mask,1))+1;
                    row=indices(index)-(column-1)*size(I_mask,1);
                    distance=((row-row_index)^2+(column-column_index)^2)^0.5;
                    if min_distance>distance && flag(index)~=1
                        min_distance=distance;
                        min_pt=index;
                    end
                end
                min_point(iter,1)=color(min_pt);
                flag(min_pt)=1;
            end
            maxfreq=mode(min_point);
            X_whole(row_index,column_index)=maxfreq;
        end
    end
    imshow(uint8(X_whole));
end