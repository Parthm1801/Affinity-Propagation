[x indices]=get_image_data_v3('s3.jpg');
 %x=importdata('x.txt');
%indices=importdata('indices.txt');
 n=length(x);
p=size(x,2);
Z=x;
m=n;
S=zeros(n);

%Computing the Similarity Matrix 
for i=1:n
    for k=1:n
        S(i,k)=-norm(x(i,:)-x(k,:)); %Similarity matrix is computed by taking negative norm of two points here, however, other methods can be used
    end
end

I=(median(median(S))*eye(n))*40; %Evaluating I using median of the off-diagonal elements of Similarity matrix
S=S+I;
S(isnan(S))=0;




% Initialize messages
N=size(S,1); 
A=zeros(N,N); %Initalizing availability matrix to zero
R=zeros(N,N); %Initalizing responsibility matrix to zero
S=S+(eps*randn(N,N))*(max(S(:))-min(S(:)))*5; % Remove degeneracies
lam=0.7; % Set damping factor


for iter=1:200 %number of iterations set to 100, can be changed depending on dataset
  
    % Compute responsibilities
    Rold=R; %Using Rold to record the previous iteration's responsibility matrix
    AS=A+S; 
    [Y,I]=max(AS,[],2); %Y stores the value of maximum value of each column of AS matrix, I stores the row index of the that value in AS
    for i=1:N 
      AS(i,I(i))=-realmax; %setting the value of maximum of AS (where i=k) equal to the IEEE minimum value
    end
    [Y2,I2]=max(AS,[],2); %second maximum value(for i=k) case
    R=S-repmat(Y,[1,N]); %repmat function makes a matrix with all N columns equal to column vector Y. Therefore, R=S-Max(A+S)
    for i=1:N 
      R(i,I(i))=S(i,I(i))-Y2(i); 
    end
    R=(1-lam)*R+lam*Rold; % Dampen responsibilities

    % Compute availabilities
    Aold=A;  %Using Aold to record the previous iteration's availibility matrix
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
idx=I(c); % Assignments, idx shows the indices of exemplars in the input data
disp(idx);



uidx=unique(idx);
for i=1:length(uidx)
    idx(uidx(i))=uidx(i);
end

disp([unique(idx),histc(idx,unique(idx))]); %displays the exemplar index and the number of elements belonging to its cluster

uidx=unique(idx);
show_image_v3(x,idx,uidx,'s3.jpg',indices)
%mkplot(x,idx);

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
       % correct_classifications=sum(color(1:100)==unique_color(1))+sum(color(101:200)==unique_color(2))+sum(color(201:300)==unique_color(3));
        correct_classifications=sum(color(1:100)==unique_color(1))+sum(color(101:200)==unique_color(2));
        if max_classifications<correct_classifications
            max_classifications=correct_classifications;
        end
    end
    error_rate=(1-max_classifications/300)*100
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
       % unique_color_all(3)=255;
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