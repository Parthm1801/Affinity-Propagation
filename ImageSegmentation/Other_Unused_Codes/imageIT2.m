% I=imread('s1.jpg');
% I=double(I);
% miny=min(min(I))*ones(size(I));
% maxy=max(max(I))*ones(size(I));
% img=(I-miny)./(maxy-miny);
% map=0*(img(:,:,1)>=0 & img(:,:,1)<0.2)+51*(img(:,:,1)>=0.2 & img(:,:,1)<0.4)+102*(img(:,:,1)>=0.4 & img(:,:,1)<0.6)+153*(img(:,:,1)>=0.6 & img(:,:,1)<0.8)+204*(img(:,:,1)>=0.8 & img(:,:,1)<1)+255*(img(:,:,1)==255);


x=get_image_data('s3.jpg',5);
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
   s(i,i)= (mins - 40*(maxs-mins));
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

disp([unique(idx),histc(idx,unique(idx))]);

uidx=unique(idx);
show_image(x,idx,uidx,'s3.jpg',5)

    
    
    function X=get_image_data(filename,interval)
    I=imread(filename);
    %I=rgb2gray(I);
    I_dash=I;
    X=[];
    for row=1:interval:floor(size(I,1)-interval)
        for column=1:interval:floor(size(I,2)-interval)
            x_unfolded=I(row:row+interval-1,column:column+interval-1,:);
            X=[X ;x_unfolded(:)'];
        end
    end
    X=double(X);
end
function show_image(X,idx,uidx,filename,interval)
%     xc=X(idx,:);
%     coX=mean(xc,2);
%     coX=(coX-min(coX))./(max(coX)-min(coX));
%   color=0*(coX(:,1)>=0 & coX(:,1)<0.2)+51*(coX(:,1)>=0.2 & coX(:,1)<0.4)+102*(coX(:,1)>=0.4 & coX(:,1)<0.6)+153*(coX(:,1)>=0.6 & coX(:,1)<0.8)+204*(coX(:,1)>=0.8 & coX(:,1)<1)+255*(coX(:,1)==255);
%   %  color=0*(coX==0)+51*(coX>0 & coX<=0.2)+102*(coX>0.2 & coX<=0.4)+153*(coX>0.4 & coX<=0.6)+204*(coX>0.6 & coX<=0.8)+255*(coX>0.8 & coX<=1);
% 
%     
%     color=idx/(sum(uidx));
%     color=color*256;
color=zeros(size(idx));
for i=1:length(uidx)
    color=color+(idx==uidx(i))*((i-1)/(length(uidx)-1));
end
color=color*255;

unique_color=unique(color);
    store=zeros(size(color,1),size(unique_color,1));
    for unique_color_index=1:size(unique_color,1)
        store(:,unique_color_index)=(unique_color(unique_color_index)==color);
    end
   [temp sorted_indices]=sort(sum(store));
   store=store(:,sorted_indices);
   store=store.*unique_color';
   color=sum(store')';
    for i=1:size(X,1)
        X(i,:)=color(i,1);
    end
    I=imread(filename);
    %I=rgb2gray(I);
    
    data_number=0;
    for row=1:interval:floor(size(I,1)-interval)
        for column=1:interval:floor(size(I,2)-interval)
            data_number=data_number+1;
            x_folded=X(data_number,:);
            I(row:row+interval-1,column:column+interval-1,:)=reshape(x_folded,interval,interval,size(I,3));
        end
    end
    I(find(125>I))=0;
    I(find(I>200))=1;
    I(find(200>=I & I>=125))=0.5;
     I=I*255;
I_mask=double(imread('s3_mask.jpg'));
    I_mask(find(125>I_mask))=0;
    I_mask(find(I_mask>200))=1;
    I_mask(find(((200>=I_mask).*I_mask)>=125))=0.5;
     I_mask=I_mask*255;
     
     %I=I/255;
     final_similarity=0;
    for unique_i=unique(I_mask)'
        indices=find(I_mask==unique_i);
        max_similarity=0;
        for unique_i_inner=unique(I_mask)'
            similarity=sum(I(indices)==unique_i_inner);
            if similarity>max_similarity
                max_similarity=similarity;
            end
        end
        final_similarity=final_similarity+max_similarity;
    end
    error_rate=(1-(final_similarity/(size(I,1)*size(I,2))))*100
    
    
    disp("show");
    imshow(I);
    
    
end

