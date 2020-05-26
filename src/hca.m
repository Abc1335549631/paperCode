function hca(input)

loadFile = input.openFile;
cluster_size = input.clusterSize;
alpha = input.alpha;
minS = input.minS;
R1 = input.R1;
lMax = input.lMax;
minSm = input.minSm;
lpH = input.lpH;
similarBnd = input.similarBnd;

%alpha=0.31736;     
%minS=59;        
%R1=0.48180;        
%lMax =46;   
%minSm=0.03319;                                                           
%similarBnd =0.1;  

data = load(loadFile);
dataNoTag = data(:,1:end-1);
dataTag = data(:,end);
xInput=dataNoTag';
Centers=[];
normalxi=[];


for i=1:size(xInput,1)
    hl=max(xInput(i,:));
    ll=min(xInput(i,:));  
    xl=(xInput(i,:)-ll)/(hl-ll);  
    normalxi=[normalxi;xl];
end
[m_data] = size(normalxi,2);
dataNorm=normalxi; 

Centers(:,1)=dataNorm(:,1);                              
Nj=1;
Ni=size(dataNorm,1);                                    


for lp=1:lpH
    cluster_cell = cell(cluster_size,1); 
    deNoise_score=zeros(1,Nj);    
    for k=1:size(dataNorm,2)                
        fireTime=zeros(Nj,1);    
            for j=1:Nj                                
                for i=1:Ni
                    fireTime(j)=fireTime(j)+(dataNorm(i,k)-Centers(i,j))*(dataNorm(i,k)-Centers(i,j));  
                end
                fireTime(j)=sqrt(fireTime(j));            

            end

        [min_fireTime,indx]=min(fireTime);                  
   
        deNoise_score(1,indx)=deNoise_score(1,indx)+1;                       

        if min_fireTime>=alpha                    
            Nj=Nj+1;                                    
            Centers=[Centers dataNorm(:,k)];         
            deNoise_score=[deNoise_score 1];           
        else

       if (k>1)||(lp>1)
            tmCenters=[];
            for j=1:Nj
                if j==indx
                     CentersJ=updateCenterJ(inf,dataNorm(:,k),Centers(:,j),Ni,inf,inf,inf); 
                else
                    CentersJ=Centers(:,j);
                end
                tmCenters=[tmCenters CentersJ];
            end   
            Centers=tmCenters;
        end
    end
  

    if Nj>1
        prn_fireTime=zeros(Nj,1);
      for j=1:Nj 
            for i=1:Ni
                prn_fireTime(j)=prn_fireTime(j)+(Centers(i,indx)-Centers(i,j))*(Centers(i,indx)-Centers(i,j));
            end 
            prn_fireTime(j)=sqrt(prn_fireTime(j));            
      end

        
        dsimiCenters=[];
        dsimi_dNscore=[];
        for j=1:Nj
          
            if (j==indx)||(prn_fireTime(j)>=similarBnd)
            
                dsimiCenters=[dsimiCenters Centers(:,j)];
                dsimi_dNscore=[dsimi_dNscore deNoise_score(:,j)];
            
            end    
        end
        
        Centers=dsimiCenters;
        deNoise_score=dsimi_dNscore;
        Nj=size(Centers,2);

    end 
    
    for cluster_num = 1:cluster_size
        if indx==cluster_num
            cluster_cell{cluster_num} = [cluster_cell{cluster_num} dataNorm(:,k)];
        end
    end   
  
end
 
% if lp>48                
%    if k==size(dataNorm,2)
%        scrThreshold=5;
%        dNcenters=[];
%    
%        for j=1:Nj
%    
%            if deNoise_score(1,j)>scrThreshold
%                dNcenters=[dNcenters Centers(:,j)];
%            end
%        end
%    
%        Centers=dNcenters;
%        Nj=size(Centers,2);
%    end
% end
    
    %dp1=[Nj k lp]  ;                                
    %dp=Centers     ;                                 
    %scr=deNoise_score;

end


                             
%c=colormap(lines(cluster_size));     
cluster_cell(cellfun(@isempty, cluster_cell))=[];
tic
cellsize=size(cluster_cell,1);
Sm = zeros(cellsize,1);
for cluster_num = 1:cellsize
     cluster = cluster_cell{cluster_num};
     a = size(cluster,2);
     if a<minS
         cluster_cell{cluster_num} = [];
         Centers(:,cluster_num)=0;
         Sm(cluster_num) = 0;
     else       
         Sm(cluster_num)=dist_4(cluster);   
     end
end
cluster_cell(cellfun(@isempty, cluster_cell))=[];     
Sm(all(Sm==0,2),:)=[];  
clusterA= cell2mat(cluster_cell');                       
try
find(ismember(dataNorm',clusterA','rows'));
catch
    addcs2
end
dataNorm=dataNorm';
dataNorm(ismember(dataNorm,clusterA','rows'),:) = [];         

Centers=Centers';
Centers(all(Centers==0,2),:)=[];    
normalxi=normalxi';
for j=1:size(Centers,1)                           
    i =0;
    c = 10;
    b = [];
for k = 1:50
    t = linspace(0,2*pi);
    ax=Centers(j,1);
    by=Centers(j,2);    
    if ax==0&&by==0
       break
    end   
    index =  j;
try
    step = exp(Sm(j))/10;  
    catch
        addcs2
end
  
    r = hypot(normalxi(:,1)-ax,normalxi(:,2)-by);                              
    T=find(r<=(R1));                           
     if T ==0 
         break 
     end 
        if i>0  

         XG=find(r<=(R1+step*((i-1)^2)*0.03));
         X1= normalxi(XG,:);

         XG1=find(r<=(R1+step*i^2*0.03));
         X2= normalxi(XG1,:);
         find(ismember(X2,X1,'rows'));      
         X2(ismember(X2,X1,'rows'),:) = [];
         idexc=find(ismember(dataNorm,X2,'rows'));
         if isempty(idexc)             
           break 
         else 
          cluster_cell{j}=[(cluster_cell{j}) dataNorm(idexc,:)']; 
          clusterA=[clusterA dataNorm(idexc,:)'];
          dataNorm(idexc,:)=[];       
          end
        end
        i=i+1;
end
    if  isempty(dataNorm)
       break
    end
end


P=[];
cellsize=size(cluster_cell,1);
for xx = 1:cellsize
    cluster = (cluster_cell{xx})';
    if ~isempty(cluster_cell{xx})
        [max_cluster]=max(cluster,[],1); 
        [min_cluster]=min(cluster,[],1);
        d1=max_cluster(1,1)-min_cluster(1,1);
        d2=max_cluster(1,2)-min_cluster(1,2) ;
        S=d1*d2;
        n=size(cluster_cell{xx},2);
        p=n/S;
        P=[P,p];                                                
    end
end
n=size(P,2);
SmL =zeros(n);
cellsize=size(cluster_cell,1);
for  i=1:n
    for j=1:n
        SmL(i,j)=abs(P(i)-P(j));                                
    end   
end  
for i = 1:n
    SmLi = SmL(i,:);
    mini = min(SmL);
    [minnum index] = find(SmL==mini);    
    dos = index(i,1);                    
   SmL(dos,dos) = inf;                                            
end
 SmL=SmL/10000;
L=zeros(n);                                                       
dc=0.25;
Dist = zeros(m_data,n);   
Centers(all(Centers==0,2),:)=[];        
for i = 1:m_data
    for j = 1:n
        Dist(i,j) = norm(normalxi(i,:)-Centers(j,:));   
    end
end
for i = 1:m_data                                                
    Disti = Dist(i,:);                              
    mini1 = min(Disti);                                           
    [minnum index] = find(Disti==mini1);           
    pos = index(1,1);         
    minL1 = minnum(1,1);      
    X = index(1,1);
    Disti(1,pos) = inf;                
    mini2 = min(Disti);                
    [minnum index] = find(Disti==mini2);   
    pos = index(1,1);
    minL2 = minnum(1,1);
    Y = index(1,1);
    if  mini2>1.5*dc
        LXY2 = 0;
    else
        LXY2 = (-1/dc)* mini2+1.5;
    end
    L(X,Y) =  L(X,Y)+LXY2;
    L(X,X) = inf;
end
  L=L+L';       
    for xx = 1:cellsize
        cluster = (cluster_cell{xx})';
        clsize=size(cluster,1);       
        lable1 = ones(clsize,1);
        lable1(:,:)=xx;
        cluster_cell{xx} = [cluster lable1];      
    end
        Number=size(normalxi,1);       
        lable1 = ones(Number,1);
        normalxi= [normalxi lable1];   
     for  xx = 1:cellsize
          cluster = (cluster_cell{xx});
          clsize=size(cluster,1);
          A=normalxi(:,1);
          B=normalxi(:,2);
          for yy=1:clsize
              x=cluster(yy,1);
              y=cluster(yy,2);
              index=find(x==A&y==B);
              normalxi(index,3)=xx;
          end
     end
Lii=[];
for i = 1:cellsize
    Li = L(i,:);
    [num indexLi] = find(Li>=lMax);
    Li(Li<lMax) = 0;
    Li(Li>=lMax) = indexLi;
    Lii = [Lii;Li];
    Lii(i,i)=i;
end  
SmLii=[];
for i = 1:cellsize
   SmLi = SmL(i,:);
   for j=1:cellsize
       if SmLi(1,j)<minSm
           SmL(i,j)=j;
       else
            SmL(i,j)=0;
       end       
   end
   SmL(i,i)=i;
end
 
   index=find((SmL- Lii)~=0);
   Lii(index)=0;
   SmL(index)=0;
   isequal( SmL,Lii)  
lable2 = 0*lable1;
normalxi= [normalxi,lable2];
C3lable = normalxi(:,end-1);                       
for i=1:cellsize
    Lii0 = Lii(i,:)';
    iv = Lii0(find(Lii0~=0),:);             
    [H_iv,~] = size(iv);
    if H_iv~=0
    miniv = min(iv);
    irm =( find(C3lable==miniv));              
    irm1 = irm(1,1);
    min0 =  normalxi(irm1,end);
    if min0==0
        min0 = i;
    end
    
    for j=1:H_iv
        k = iv(j,1);
        irC3lable = find(C3lable==k);
        [h_irC3lable,~]=size(irC3lable);
        for tt=1:h_irC3lable
            tt0=irC3lable(tt,1);
            normalxi(tt0,end) = min0;
        end
    end
    end
end

[m_C3,L_C3] = size(normalxi);
C4 = sortrows(normalxi,L_C3);

C5 = sortrows(normalxi,L_C3-2);  %? not use;
class = C4(:,end);
c = max(class);
x1 = C4(:,1:end-2);
x1=[x1 class];
cmap=colormap; 
 
for i=1:m_C3
    if x1(i,3)~=-1
   ic=int8((x1(i,3)*100.)/(c*10));

    plot(x1(i,1),x1(i,2),'o','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));   
    hold on;
    end
    if x1(i,3)==-1
        plot(x1(i,1),x1(i,2),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
         hold on;
    end
end
toc
  [FMeasure2,Accuracy2] = Fmeasure(dataTag,normalxi(:,end))  

  
end


