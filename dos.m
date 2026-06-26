n=size(siteseta,2);
%tt=9;
ifnall=[];
infeall=[];
pos=[];
site=[];
oder=[];
ms=0;
cof=2;
xi=1:0.5:tt;

for i=1:n
    [siter,siter_aa]=oringnal_data(siteseta(i));
    site1=siter(1,1:tt);%interp1(1:tt,siter(1,1:tt), xi, 'linear');%
    site2=siter(2,1:tt);%interp1(1:tt,siter(2,1:tt), xi, 'linear');%
    site3=siter(3,1:tt);%interp1(1:tt,siter(3,1:tt), xi, 'linear');%
    site4=siter(4,1:tt);
    siter_aa=siter_aa(:,tt)';
    fo=emd(site1);
    ifnall2=[];
    infeall2=[];
    if isempty(fo)==0
        [hs,f,t,imfinsf,imfinse]=hht(fo);
        marginal_spectrum = sum(hs, 2); 
        ms=ms+marginal_spectrum;
        imfinse(:,imfinsf(tt,:)<0)=[];
        imfinsf(:,imfinsf(tt,:)<0)=[];
        [s,ss]=size(imfinsf);
        [aa,bb]=max(imfinsf(tt,:));
        if ss>1 && aa>cof
            ss=ss-1;
            imfinse(:,bb)=[];
            imfinsf(:,bb)=[];        
        end
        ifnall2=imfinsf;
        infeall2=imfinse;
        oder=[oder,siter_aa(1)*ones([1,ss])];
        %site=[site,siteseta(i)];
    end
    if site2(tt)>0.3*site1(tt)
        fo=emd(site2);
        %[s,ss]=size(fo);
        if isempty(fo)==0
            [hs,f,t,imfinsf,imfinse]=hht(fo);
            marginal_spectrum = sum(hs, 2); 
            ms=ms+marginal_spectrum;
            imfinse(:,imfinsf(tt,:)<0)=[];
            imfinsf(:,imfinsf(tt,:)<0)=[];
            [s,ss]=size(imfinsf);
            [aa,bb]=max(imfinsf(tt,:));
            if ss>1 && aa>cof
               ss=ss-1;
               imfinse(:,bb)=[];
               imfinsf(:,bb)=[];        
            end
            infeall2=[infeall2,imfinse];
            ifnall2=[ifnall2,imfinsf];
            oder=[oder,siter_aa(2)*ones([1,ss])];
            %site=[site,siteseta(i)];
        end
    end
    if site2(tt)>0.3*site1(tt)&& site3(tt)>0.3*site2(tt)
        fo=emd(site3);
        if isempty(fo)==0
            [hs,f,t,imfinsf,imfinse]=hht(fo);
            marginal_spectrum = sum(hs, 2); 
            ms=ms+marginal_spectrum;
            imfinse(:,imfinsf(tt,:)<0)=[];
            imfinsf(:,imfinsf(tt,:)<0)=[];
            [s,ss]=size(imfinsf);
            [aa,bb]=max(imfinsf(tt,:));
            if ss>1 && aa>cof
               ss=ss-1;
               imfinse(:,bb)=[];
               imfinsf(:,bb)=[];        
            end
            infeall2=[infeall2,imfinse];
            ifnall2=[ifnall2,imfinsf];
            oder=[oder,siter_aa(3)*ones([1,ss])];
            %site=[site,siteseta(i)];
        end
    end
    
    if site2(tt)>1*site1(tt)&& site4(tt)>0*site3(tt)
        fo=emd(site4);
        if isempty(fo)==0
            [hs,f,t,imfinsf,imfinse]=hht(fo);
            marginal_spectrum = sum(hs, 2); 
            ms=ms+marginal_spectrum;
            imfinse(:,imfinsf(tt,:)<0)=[];
            imfinsf(:,imfinsf(tt,:)<0)=[];
            [s,ss]=size(imfinsf);
            [aa,bb]=max(imfinsf(tt,:));
            if ss>1 && aa>cof
               ss=ss-1;
               imfinse(:,bb)=[];
               imfinsf(:,bb)=[];        
            end
            infeall2=[infeall2,imfinse];
            ifnall2=[ifnall2,imfinsf];
            oder=[oder,siter_aa(4)*ones([1,ss])];
            %site=[site,siteseta(i)];
        end
    end
    
    %if isempty(ifnall2)==0 && sum(ifnall2(tt,:)>a0 &ifnall2(tt,:)<b)>0
        %pos=[i,pos];
     %   pos=[siteseta(i),pos];
    %end
    %ifnall2(ifnall2>2.5)=[];
    
    [s,ss]=size(ifnall2);
    s=siteseta(i)*ones([1,ss]);
    site=[site,s];
    
    %(ifnall2(tt,:))
    infeall=[infeall,infeall2];
    ifnall=[ifnall,ifnall2];
end
edges = linspace(0, 3, 100);
pt=histogram(ifnall(tt,:),edges);
xlim([0,3])
xlabel('f (frequency)','FontSize', 15)
ylabel('number of modes','FontSize', 15)
unique(pos);
%title('2020-09','FontSize',15)