%[x indices]=getfile('wolf.jpg',20);
x=importdata('xWolfCrisp.txt');
indices=importdata('WolfCrispIndex.txt');
% n=length(x);
% p=size(x,2);
% Z=x;
% m=n;
% 
% 
% %Computing the Similarity Matrix 
% for i=1:n
%     for k=1:n
%         S(i,k)=-norm(x(i,:)-x(k,:)); %Similarity matrix is computed by taking negative norm of two points here, however, other methods can be used
%     end
% end
% 
% I=(median(median(S))*eye(n))*500; %Evaluating I using median of the off-diagonal elements of Similarity matrix
% S=S+I;
% S(isnan(S))=0;
% 
% 
% 
% 
% % Initialize messages
% N=size(S,1); 
% A=zeros(N,N); %Initalizing availability matrix to zero
% R=zeros(N,N); %Initalizing responsibility matrix to zero
% S=S+(eps*randn(N,N))*(max(S(:))-min(S(:)))*5; % Remove degeneracies
% lam=0.9; % Set damping factor
% 
% 
% for iter=1:200 %number of iterations set to 100, can be changed depending on dataset
%   
%     % Compute responsibilities
%     Rold=R; %Using Rold to record the previous iteration's responsibility matrix
%     AS=A+S; 
%     [Y,I]=max(AS,[],2); %Y stores the value of maximum value of each column of AS matrix, I stores the row index of the that value in AS
%     for i=1:N 
%       AS(i,I(i))=-realmax; %setting the value of maximum of AS (where i=k) equal to the IEEE minimum value
%     end
%     [Y2,I2]=max(AS,[],2); %second maximum value(for i=k) case
%     R=S-repmat(Y,[1,N]); %repmat function makes a matrix with all N columns equal to column vector Y. Therefore, R=S-Max(A+S)
%     for i=1:N 
%       R(i,I(i))=S(i,I(i))-Y2(i); 
%     end
%     R=(1-lam)*R+lam*Rold; % Dampen responsibilities
% 
%     % Compute availabilities
%     Aold=A;  %Using Aold to record the previous iteration's availibility matrix
%     Rp=max(R,0); 
%     for k=1:N 
%       Rp(k,k)=R(k,k); 
%     end
%     A=repmat(sum(Rp,1),[N,1])-Rp;
%     dA=diag(A); 
%     A=min(A,0);
%     for k=1:N 
%       A(k,k)=dA(k);
%     end
%     A=(1-lam)*A+lam*Aold; % Dampen availabilities
%     
% end


r=load('respCrispWolf.mat');
R=r.R;
a=load('avCrispWolf.mat');
A=a.A;
s=load('simCrispWolf.mat');
S=s.S;


E=R+A; % Pseudomarginals
I=find(diag(E)>0); 
K=length(I); % Indices of exemplars
[tmp,c]=max(S(:,I),[],2); 
idx=I(c); % Assignments, idx shows the indices of exemplars in the input data




uidx=unique(idx);
for i=1:length(uidx)
    idx(uidx(i))=uidx(i);
end

disp([unique(idx),histc(idx,unique(idx))]); %displays the exemplar index and the number of elements belonging to its cluster

uidx=unique(idx);
create_image_2(x,idx,uidx,indices,'wolf.jpg');




 
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
    color(find(color==1))=0;
    color(find(color==0.5))=1;
    color(find(color==13))=0.5;
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

