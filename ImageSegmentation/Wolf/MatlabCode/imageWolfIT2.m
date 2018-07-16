
%[x indices]=getfile('wolf.jpg',20);
x=importdata('x_Wolf.txt');
indices=importdata('wolfindices.txt');
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
beta1=1;
beta2=10;
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
disbar(isnan(disbar))=0;

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
   s(i,i)= (mins - 50*(maxs-mins));
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
create_image_2(x,idx,uidx,indices,'wolf.jpg')




 
function [X all_indices]=getfile(filename,n)
I=imread(filename);
   
       all_indices=ceil(rand(1,ceil(size(I,1)*size(I,2)/n))*(size(I,1)*size(I,2)))';
       I1=I(:,:,1);
       I2=I(:,:,2);
       I3=I(:,:,3);
       I=rgb2gray(I);
       I1=imgaussfilt(I,0.5);%standard deviation set to 0.5 gaussian blur
       I2=medfilt2(I,[3 3]);%m*n box is 3*3; median iltering
       I3_temp=entropyfilt(I,true((9)));%entropy filtering
       I3=(I3_temp-min(I3_temp))./(max(I3_temp)-min(I3_temp));
       X=[I1(all_indices) I2(all_indices) I3(all_indices)];
       X=double(X);
        return;
    end
      
function create_image_2(X,idx,uidx,indices,filename)
    Img=imread(filename);
    Img=rgb2gray(Img);
    color=zeros(size(idx));
    for i=1:length(uidx)
        color=color+(idx==uidx(i))*((i-1)/(length(uidx)-1));
    end
        color(find(color==0))=13;
    color(find(color==0.5))=0;
    color(find(color==1))=0.5;
    color(find(color==13))=1;
    color=color*255;  
    X_whole=zeros(size(Img));    
    for row_index=1:size(Img,1)
        for column_index=1:size(Img,2)
            min_distance=inf;
            min_point=zeros(5,1);
            flag=zeros(size(indices,1),1);
            for iter=1:5
                for index=1:size(indices,1)
                    column=floor((indices(index)-0.1)/size(Img,1))+1;
                    row=indices(index)-(column-1)*size(Img,1);
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
    X_whole=uint8(X_whole);
    
    
    I_mask=(imread('wolf_mask.jpg'));
    I_mask=double(rgb2gray(I_mask));
    I_mask(find(125>I_mask))=0;
    I_mask(find(I_mask>200))=1;
    I_mask(find(((200>=I_mask).*I_mask)>=125))=0.5;
    I_mask=uint8(I_mask*255);
    max_classifications=0;
    unique_color_all=unique(color);
    negative_indices=0;
    indices_list=cell(size(unique_color_all));
    unique_color=unique(X_whole)';
    for unique_color_index=1:size(unique_color,2)
        indices=find(X_whole==unique_color(unique_color_index));
        indices_list(unique_color_index)={indices};
    end
    max_classifications=0;
    confusion_matrix_2=zeros(size(unique_color,1),size(unique_color,1));
    if sum(unique_color_all>=0)
        unique_color_all=zeros(3,1);
        unique_color_all(1)=0;
        unique_color_all(2)=128;
        unique_color_all(3)=255;
        
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
   confusion_matrix=zeros(3); 
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

type=['s' 'o' 'x'];
for j=1:3
X(:,j)=nrmlz(X(:,j));
end
plotcube([1 1 1],[0 0 0],0);
hold on
for i=1:length(uidx)
ind=find(idx==uidx(i));
scatter3(X(ind,1),X(ind,2),X(ind,3),100,color(ind),type(i));
end
xlabel('Entropy','FontSize',14,'FontWeight','bold');
ylabel('Median Filtered Intensity','FontSize',14,'FontWeight','bold');
zlabel('Gaussian Blur','FontSize',14,'FontWeight','bold');
set(gca,'FontName','Times New Roman');
end
 
 
 function Nx=nrmlz(x)
    Nx=(x-min(x))./(max(x)-min(x));
 end
 
 
