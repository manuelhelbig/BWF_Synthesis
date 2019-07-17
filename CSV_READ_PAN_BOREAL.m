%%% read csv files for EC sites %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all
clc
Lv         = 2.5e6;                 % heat of vaporization [J kg-1]
% read site data
%filename = '/Users/manuelhelbig/Desktop/Mac_PDF_Boreal/PFT_SITES.csv';
cd('E:\Mac_PDF_Boreal\');
filename = 'E:\Mac_PDF_Boreal\PFT_SITES_4.csv';
delimiter = ',';
startRow = 2;

formatSpec = '%s%s%f%f%f%f%f%f%[^\n\r]';

fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
% Allocate imported array to column variable names
SITE = dataArray{:, 1};
PFT = dataArray{:, 2};
DIST = dataArray{:, 3};
LAT = dataArray{:, 4};
LONG = dataArray{:, 5};
TMP = dataArray{:, 6};
PRE = dataArray{:, 7};
PET = dataArray{:, 8};

clearvars filename delimiter startRow formatSpec fileID dataArray ans;
%%
%cd('/Users/manuelhelbig/Desktop/Mac_PDF_Boreal');
%path='/Users/manuelhelbig/Desktop/Mac_PDF_Boreal/';
path='E:\Mac_PDF_Boreal\';

finfo = dir([path '*.csv']);
n=0;
lng=90;
for k = 1:length(finfo);
    if strcmp(finfo(k).name(3),'-') && (length(finfo(k).name)==10 | length(finfo(k).name)==11);
        n=n+1;
        finfo(k).name
        SITE_NAME{n}=finfo(k).name(1:6);
    end
end
%%
DT_ALL_YR=datenum(1993,1,1,0,30,0):1/48:datenum(2018,12,31,23,30,0);
DT_ALL_YR=floor(DT_ALL_YR.*48)./48;
COVERAGE=NaN(lng,length(DT_ALL_YR));
MM_LNGT=[31 28 31 30 31 30 31 31 30 31 30 31];
%%
for k = 1:lng; 
    SITE_NAME{k}
    data=importdata([SITE_NAME{k} '.csv']);
   
    if k == 1;
        LABEL=strsplit(data.textdata{1,1},',');
    end
    i1=find(strcmp(LABEL,'VPD'));
    i2=find(strcmp(LABEL,'LE'));
    i3=find(strcmp(LABEL,'Rnet'));
    i4=find(strcmp(LABEL,'DOY'));
    i5=find(strcmp(LABEL,'LE_GF'));
    i6=find(strcmp(LABEL,'YEAR'));
    i7=find(strcmp(LABEL,'HOUR'));
    i8=find(strcmp(LABEL,'MIN'));
    iTa=find(strcmp(LABEL,'TA'));
    TA_SITES{k}=data.data(:,iTa);
    %{
    VPD_LE{k}(:,1)=data.data(:,i1);
    VPD_LE{k}(:,2)=data.data(:,i2);
    
    
    if strcmp(SITE_NAME{k},'SE-RO2') | strcmp(SITE_NAME{k},'SE-RO3')
        data.data(:,i5)=data.data(:,i5).*30.*60.*(18.01./1000)./1000;
    end
    RN=data.data(:,i3);
    RN(RN==0)=NaN;
    VPD_LE{k}(:,3)=RN;
    clear RN
    VPD_LE{k}(:,4)=data.data(:,i4);
    VPD_LE{k}(:,5)=data.data(:,i5);
    VPD_LE{k}(:,6)=data.data(:,i7);
    VPD_LE{k}(:,7)=data.data(:,i8);
    
    i9=find(strcmp(LABEL,'H'));
    if length(data.data(1,:))<21
        VPD_LE{k}(:,8)=NaN;
    else
        VPD_LE{k}(:,8)=data.data(:,i9);
    end
    
    % check for periods of high Rnet and low LE
    if length(data.data(1,:))<21
        GS_DERIV(k)=0;
        VPD_LE{k}(:,9)=NaN;
        VPD_LE{k}(:,10)=NaN;
    elseif sum(~isnan(data.data(:,find(strcmp(LABEL,'USTAR')))))>0 & sum(~isnan(data.data(:,find(strcmp(LABEL,'WS')))))>0 ...
            & sum(~isnan(data.data(:,find(strcmp(LABEL,'TA')))))>0 & sum(~isnan(data.data(:,find(strcmp(LABEL,'VPD')))))>0 ...
            & sum(~isnan(data.data(:,find(strcmp(LABEL,'H')))))>0 & sum(~isnan(data.data(:,find(strcmp(LABEL,'TA')))))>0
        GS_DERIV(k)=1;
        gaero = ga(data.data(:,find(strcmp(LABEL,'USTAR'))),data.data(:,find(strcmp(LABEL,'WS'))));
        for s=1:3;
            gaero=filter_papale(gaero,isnan(gaero),4.5,500,VPD_LE{k}(:,1),600);
        end
        close all
        gaero(gaero<-1 | gaero==0)=NaN;
        FILT1=gaero<prctile(gaero,0.5);
        FILT2=gaero>prctile(gaero,99.5);
        gaero(FILT1==1)=NaN;
        gaero(FILT2==1)=NaN;
        clear FILT1 FILT2

        [gsurf] = gs(VPD_LE{k}(:,8),VPD_LE{k}(:,2),...
            VPD_LE{k}(:,1),data.data(:,find(strcmp(LABEL,'TA'))),...
            gaero,0,0,1,NaN);
        for s=1:3;
            gsurf=filter_papale(gsurf,isnan(gsurf),4.5,500,VPD_LE{k}(:,1),600);
        end
        close all
        gsurf(gsurf<-1 | gsurf==0)=NaN;
        FILT1=gsurf<prctile(gsurf,0.5);
        FILT2=gsurf>prctile(gsurf,99.5);
        gsurf(FILT1==1)=NaN;
        gsurf(FILT2==1)=NaN;
        clear FILT1 FILT2
        
        VPD_LE{k}(:,9)=gsurf;
        VPD_LE{k}(:,10)=gaero;
        clear gs ga
    else
        GS_DERIV(k)=0;
        VPD_LE{k}(:,9)=NaN;
        VPD_LE{k}(:,10)=NaN;
    end
    i10=find(strcmp(LABEL,'SWin'));
    i11=find(strcmp(LABEL,'G'));
    if i11<=length(data.data(1,:));
        VPD_LE{k}(:,11)=data.data(:,i11);
        
    else
        VPD_LE{k}(:,11)=NaN;
    end
    VPD_LE{k}(:,12)=data.data(:,i10);
    %}
    YR{k}=data.data(:,i6);
    %{
    i = find(strcmp(SITE,SITE_NAME{k}));
    if strcmp(SITE_NAME{k},'RU-YLF');
        PFT_NUM{k}='DNF';
        DIST_NUM(k)=0;
        PRE_NUM(k)=PRE(find(strcmp(SITE,'RU-SKP')));
        PET_NUM(k)=NaN;
        LAT_NUM(k)=LAT(find(strcmp(SITE,'RU-SKP')));
        LON_NUM(k)=LONG(find(strcmp(SITE,'RU-SKP')));
    else
        PFT_NUM{k}=PFT{i};
        DIST_NUM(k)=DIST(i);
        PRE_NUM(k)=PRE(i);
        PET_NUM(k)=PET(i);
        LAT_NUM(k)=LAT(i);
        LON_NUM(k)=LONG(i);
    end
    clear i
    
    % get monthly data
    YR=data.data(:,i6);
    DOY=floor(data.data(:,i4));
    DT_TEMP=doy2date(DOY,YR)+VPD_LE{k}(:,6)./24+(VPD_LE{k}(:,7)./60)./24;
    DT_TEMP=floor(DT_TEMP.*48)./48;
    sel_COV=find(ismember(DT_ALL_YR,DT_TEMP));
    %sel_COV_i=find(ismember(DT_TEMP,DT_ALL_YR));
    if strcmp(SITE_NAME{k},'CA-NSF');
        sel_NSF=vertcat(1,find(DT_TEMP(2:end)-DT_TEMP(1:end-1)~=0)+1);
        COVERAGE(k,sel_COV)=VPD_LE{k}(sel_NSF,2);
    elseif strcmp(SITE_NAME{k},'RU-YLF');
        for j = 1:length(DT_TEMP);
            sel_j=find(DT_ALL_YR==DT_TEMP(j));
            if ~isempty(sel_j);
                COVERAGE(k,sel_j)=VPD_LE{k}(j,2);
            end
            clear sel_j
        end
    else
        COVERAGE(k,sel_COV)=VPD_LE{k}(:,2);
    end
    clear sel_COV sel_COV_i DT_TEMP sel_NSF
    %subplot(2,1,1)
    %plot(doy2date(data.data(:,i4),YR),VPD_LE{k}(:,2),'.','MarkerSize',4)
    %hold on
    MM=str2num(datestr(doy2date(DOY,YR),'mm'));
    uYR=unique(YR);
    HOUR=data.data(:,i7);
    w=0;
    for r=1:length(unique(YR));
        % choose entire year
        ind=find(YR==uYR(r));
        % choose growing season
        ind1=find(YR==uYR(r) & (DOY>=121 & DOY<=273));
        MM_IND=MM(ind);
        uMM=unique(MM_IND);
        LE_MM_GF=VPD_LE{k}(ind,5);
        LE_MM_GF1=VPD_LE{k}(ind1,5);
        LE_MM=VPD_LE{k}(ind,2);
        for l=1:length(uMM);
            w=w+1;
            i2=find(MM_IND==uMM(l));
            MEAN_LE(w)=nanmean(LE_MM_GF(i2));
            
        
            if sum(~isnan(LE_MM_GF(i2)))/(48.*MM_LNGT(uMM(l)))>=0.9 & sum(~isnan(LE_MM(i2)))/(48.*MM_LNGT(uMM(l)))>=0.3
                   
                if prctile(VPD_LE{k}(:,5),95)>5;
                    
                   LE_MNT(w)=nanmean(LE_MM_GF(i2).*30.*60./Lv).*48.*MM_LNGT(uMM(l));
                else
                   LE_MNT(w)=nanmean(LE_MM_GF(i2)).*48.*MM_LNGT(uMM(l));
                end
            else
                LE_MNT(w)=NaN;
            end
            VPD_MM=VPD_LE{k}(ind,1);
            if sum(~isnan(VPD_MM(i2)))/(48.*MM_LNGT(l))>=0.8
                i = find(YR==uYR(r) & MM==uMM(l) & (HOUR>=12 & HOUR<=15));
                VPD_MNT(w)=nanmean(VPD_MM(i2));
                VPD_MNT_MID(w)=nanmean(VPD_LE{k}(i,1));
            else
                VPD_MNT(w)=NaN;
                VPD_MNT_MID(w)=NaN;
            end
            clear i
            Rn_MM=VPD_LE{k}(ind,3);
            Rn_MNT(w)=nanmean(Rn_MM(i2));
            SW_MM=data.data(ind,find(strcmp(LABEL,'SWin')));
            SW_MNT(w)=nanmean(SW_MM(i2));
            TA_MM=data.data(ind,find(strcmp(LABEL,'TA')));
            TA_MNT(w)=nanmean(TA_MM(i2));
            YY_MNT(w)=uYR(r);
            MNT(w) = uMM(l);
            PFT_MNT(w)=double(strcmp(PFT_NUM{k},'WET'));
        end
        clear TA_MM SW_MM LE_MM LE_GF_MM Rn_MM
    end  
    
%     plot(VPD_MNT_MID,LE_MNT,'o');
%     title(SITE_NAME{k});
%     pause()
%     close all
    YR=data.data(:,i6);
    DOY=floor(data.data(:,i4));
    uYR=unique(YR);
    HOUR=data.data(:,i7);
    for r=1:length(unique(YR));
        ind=find(YR==uYR(r) & (DOY>=121 & DOY<=273));
        MEAN_LE(r)=nanmean((VPD_LE{k}(ind,5)));
        
        GAP = 0;
        maxGAP = 0;
        for s = 1:length(ind);
            if isnan(VPD_LE{k}(ind(s),5))
                GAP = GAP+1;
            else
                GAP = 0;
            end
        
            if GAP>maxGAP;
                maxGAP=GAP;
            end
        end
        clear GAP
        
        if sum(~isnan(VPD_LE{k}(ind,5)))/(48.*153)>=0.9 & sum(~isnan(VPD_LE{k}(ind,2)))/(48.*153)>=0.3 & maxGAP<21*48
            
            if nanmean((VPD_LE{k}(ind,5)))>5;
                LE_GS(r)=nanmean(VPD_LE{k}(ind,5).*30.*60./Lv).*48.*153;
            else
                LE_GS(r)=nanmean(VPD_LE{k}(ind,5)).*48.*153;
            end
        else
            LE_GS(r)=NaN;
        end
        clear maxGAP
        
        if sum(~isnan(VPD_LE{k}(ind,1)))/(48.*153)>=0.8
            i = find(YR==uYR(r) & (DOY>=121 & DOY<=273) & (HOUR>=12 & HOUR<=15));
            VPD_GS(r)=nanmean(VPD_LE{k}(ind,1));
            VPD_MID(r)=nanmean(VPD_LE{k}(i,1));
        else
            VPD_GS(r)=NaN;
            VPD_MID(r)=NaN;
        end
        clear i
        Rn_GS(r)=nanmean(VPD_LE{k}(ind,3));
        SWin_GS(r)=nanmean(data.data(ind,find(strcmp(LABEL,'SWin'))));
        TA_GS(r)=nanmean(data.data(ind,find(strcmp(LABEL,'TA'))));
        YY(r)=uYR(r);
        PFT_ANN(r)=double(strcmp(PFT_NUM{k},'WET'));
        PRE_ANN(r)=PRE_NUM(k);
        LAT_ANN(r)=LAT_NUM(k);
        LON_ANN(r)=LON_NUM(k);
        clear ind
        
        ind=find(YR==uYR(r));
                
        GAP = 0;
        maxGAP = 0;
        for s = 1:length(ind);
            if isnan(VPD_LE{k}(ind(s),5))
                GAP = GAP+1;
            else
                GAP = 0;
            end
        
            if GAP>maxGAP;
                maxGAP=GAP;
            end
        end
        clear GAP
        %maxGAP
        if sum(~isnan(VPD_LE{k}(ind,5)))/(48.*365)>=0.8 & sum(~isnan(VPD_LE{k}(ind,2)))/(48.*365)>=0.25 & maxGAP<21*48
            if nanmean((VPD_LE{k}(ind,5)))>5;
                LE_ANN(r)=nanmean(VPD_LE{k}(ind,5).*30.*60./Lv).*48.*365;
            else
                LE_ANN(r)=nanmean(VPD_LE{k}(ind,5)).*48.*365;
            end
        else
            LE_ANN(r)=NaN;
        end
        clear ind
    end
    
    clear maxGAP
%     MEAN_LE
%     LE_GS
%         plot(data.data(:,i5),'.');
%         hold on
%         plot(data.data(:,i2),'.');
%         pause()
%         close all
%         
%         plot(data.data(:,i1),'.');
%         pause()
%         close all
    
%}    
    
%     MAT=horzcat(YY_MNT',MNT',LE_MNT',VPD_MNT',VPD_MNT_MID',Rn_MNT',SW_MNT',TA_MNT',PFT_MNT');
%     eval(['csvwrite(''MNT_' SITE_NAME{k} '.csv'',MAT);']);
%     clear MAT
     clear YY_MNT MNT LE_MNT VPD_MNT VPD_MNT_MID Rn_MNT SW_MNT TA_MNT PFT_MNT
%     MAT=horzcat(YY',LE_ANN',LE_GS',VPD_GS',Rn_GS',SWin_GS',TA_GS',PFT_ANN',PRE_ANN',LAT_ANN',LON_ANN',VPD_MID');
%     eval(['csvwrite(''ANN_' SITE_NAME{k} '.csv'',MAT);']);
    clear MAT YY LE_ANN LE_GS VPD_GS Rn_GS SWin_GS TA_GS PFT_ANN MEAN_LE PRE_ANN LAT_ANN LON_ANN PET_ANN VPD_MID  
end
%%
figure,
subplot(2,1,1);
plot(DT_ALL_YR,COVERAGE,'.','MarkerSize',3);
xlim([DT_ALL_YR(1) DT_ALL_YR(end)]);
ylim([-50 650]);
ylabel('Latent heat flux (LE) [W m^{-2}]');
subplot(2,1,2);
plot(DT_ALL_YR,sum(~isnan(COVERAGE)),'.','MarkerSize',3);
xlim([DT_ALL_YR(1) DT_ALL_YR(end)]);
datetickzoom('x','yyyy');
ylabel('#sites with LE');
%%
VPD_BINS=NaN(lng,40);
LE_BINS=NaN(lng,40);

VPD_BINS_SIT=NaN(lng,40);
LE_BINS_SIT=NaN(lng,40);

% see Oren 1999 (accounts for decrease in stomatal conductance with VPD)
VPD_RESP=@(par,VPD) ((-par(1).*log(VPD)+par(2)).*VPD);
    opts = statset('nlinfit');
    opts.RobustWgtFun = 'bisquare';
    
for k=1:lng;
    ind = find((VPD_LE{k}(:,4)<121 | VPD_LE{k}(:,4)>=274));% & (VPD_LE{k}(:,6)<=12 | VPD_LE{k}(:,6)>=16));   
    %ind = find((VPD_LE{k}(:,4)<152 | VPD_LE{k}(:,4)>=196));
    %ind = find((VPD_LE{k}(:,4)<196 | VPD_LE{k}(:,4)>=244));
    VPD_GS=VPD_LE{k}(:,1);
    VPD_GS(ind)=NaN;
    Rnet_GS=VPD_LE{k}(:,3);
    Rnet_GS(ind)=NaN;
    LE_GS=VPD_LE{k}(:,2).*30.*60./Lv.*2;
    LE_GS(ind)=NaN;
    clear ind
    ind=find(~isnan(VPD_GS)&~isnan(LE_GS));
    %ind=find(~isnan(VPD_GS)&~isnan(LE_GS)&~isnan(Rnet_GS));
    x=VPD_GS(ind);
    y=LE_GS(ind);
    
    z=Rnet_GS(ind);
    
    clear ind
    
    [par]=polyfit(z,y,1);
%     figure,
%     plot(z,y,'.');
%     hold on
%     plot(z,par(1).*z+par(2),'-');
%     title(SITE_NAME{k});
%     pause()
%     close all
    pr(:,k) = nlinfit(x,y,VPD_RESP,[0.05 0.1],opts); 
   %y=y-par(1).*z+par(2);
   if PFT_NUMERIC(k)==1;
        plot(min(x):0.1:max(x),VPD_RESP(pr(:,k),min(x):0.1:max(x)),'r','LineWidth',2);
   else
       plot(min(x):0.1:max(x),VPD_RESP(pr(:,k),min(x):0.1:max(x)),'b','LineWidth',2);
   end
   hold on
   xlim([0 5]);
   ylim([0 0.4]);
   axis square
   title(SITE_NAME{k})
   %pause(3)
   %close all
    PRC_VPD_SIT=prctile(x,[2.5:2.5:97.5]);
    PRC_VPD=0.1:0.075:3;
    for l =1:40;
        if l == 1;
            ind=find(x<=PRC_VPD(l));
            ind_SIT=find(x<=PRC_VPD_SIT(l));
        elseif l<=39;
            ind=find(x>PRC_VPD(l-1) & x<=PRC_VPD(l));
            ind_SIT=find(x>PRC_VPD_SIT(l-1) & x<=PRC_VPD_SIT(l));
        else
            ind=find(x>PRC_VPD(l-1));
            ind_SIT=find(x>PRC_VPD_SIT(l-1));
        end
        if length(ind)<3 | isempty(ind)
            VPD_BINS(k,l)=NaN;
            LE_BINS(k,l)=NaN;
        else
            VPD_BINS(k,l)=nanmean(x(ind));
            LE_BINS(k,l)=nanmean(y(ind));
        end
        VPD_BINS_SIT(k,l)=nanmean(x(ind_SIT));
        LE_BINS_SIT(k,l)=nanmean(y(ind_SIT));
    end
    clear PRC_VPD x y VPD_GS LE_GS ind ind_SIT
end

for k=1:90;
    if strcmp(PFT_NUM{k},'WET');
        PFT_NUMERIC(k)=1;
    else
        PFT_NUMERIC(k)=0;
    end
end
% for k = 1:84;
%     LE_BINS(k,:)=LE_BINS(k,:)./LE_BINS(k,13);
% end

%% derive surface conductance sensitivity to VPD
VPD_BINS=NaN(lng,40);
GS_BINS=NaN(lng,40);

VPD_BINS_SIT=NaN(lng,40);
GS_BINS_SIT=NaN(lng,40);
for k=1:lng;
    if GS_DERIV(k)==0;
        VPD_BINS(k,:)=NaN;
        GS_BINS(k,:)=NaN;
    else
        % filter data (growing season, LE > 50 Wm-2, between 6am and 9pm,
        % and LE+H>100 Wm-2)
        DOY_start=182;
        DOY_end=212;
        ind = find((VPD_LE{k}(:,4)<DOY_start | VPD_LE{k}(:,4)>=DOY_end) | VPD_LE{k}(:,2)<50 ...
            | (VPD_LE{k}(:,6)<9 | VPD_LE{k}(:,6)>21) | (VPD_LE{k}(:,8)+VPD_LE{k}(:,2))<100);
        VPD_GS=VPD_LE{k}(:,1);
        VPD_GS(ind)=NaN;
        Rnet_GS=VPD_LE{k}(:,3);
        Rnet_GS(ind)=NaN;
        GS_GS=VPD_LE{k}(:,9);
        GS_GS(ind)=NaN;
        clear ind
        % only take positive gs
        i=find(~isnan(VPD_GS)&~isnan(GS_GS)&GS_GS>0);
        x=VPD_GS(i);
        y=GS_GS(i);
    
        clear ind
    
        PRC_VPD_SIT=prctile(x,[2.5:2.5:97.5]);
        PRC_VPD=0.1:0.075:3;
        for l =1:40;
            if l == 1;
                ind=find(x<=PRC_VPD(l));
                ind_SIT=find(x<=PRC_VPD_SIT(l));
            elseif l<=39;
                ind=find(x>PRC_VPD(l-1) & x<=PRC_VPD(l));
                ind_SIT=find(x>PRC_VPD_SIT(l-1) & x<=PRC_VPD_SIT(l));
            else
                ind=find(x>PRC_VPD(l-1));
                ind_SIT=find(x>PRC_VPD_SIT(l-1));
            end
            if length(ind)<3 | isempty(i)
                VPD_BINS(k,l)=NaN;
            	GS_BINS(k,l)=NaN;
                VPD_BINS_UP(k,l)=NaN;
                GS_BINS_UP(k,l)=NaN;
            else
                VPD_BINS(k,l)=nanmean(x(ind));
                GS_BINS(k,l)=nanmean(y(ind));
                x1=x(ind);
                y1=y(ind);
                VPD_BINS_UP(k,l)=nanmean(x1>prctile(x1,80));
                GS_BINS_UP(k,l)=nanmean(y1>prctile(y1,80));
            end
            VPD_BINS_SIT(k,l)=nanmean(x(ind_SIT));
            GS_BINS_SIT(k,l)=nanmean(y(ind_SIT));
        end
    clear x y VPD_GS GS_GS ind ind_SIT i
    end
end
i=find(strcmp(PFT_NUM,'WET'));
i1=find(~strcmp(PFT_NUM,'WET'));
gaero_WET=[];
gaero_FOR=[];
for k=1:lng;
    
    gaero=VPD_LE{k}(:,10);
    ind = find((VPD_LE{k}(:,4)>=DOY_start & VPD_LE{k}(:,4)<DOY_end) & (VPD_LE{k}(:,6)>=9 & VPD_LE{k}(:,6)<=21));
    gaero=gaero(ind);
    if ismember(k,i);
        gaero_WET=vertcat(gaero_WET,gaero(~isnan(gaero)));
    else
        gaero_FOR=vertcat(gaero_FOR,gaero(~isnan(gaero)));
    end
    clear gaero
end
%%
gaero=NaN(length(gaero_FOR),2);
gaero(:,1)=gaero_FOR;
gaero(1:length(gaero_WET),2)=gaero_WET;
figure,
iosr.statistics.boxPlot(gaero,'showOutliers',false)
%% calculate if GS per VPD is significantly different
i=find(strcmp(PFT_NUM,'WET'));
i1=find(~strcmp(PFT_NUM,'WET'));
map_CB=flipud(brewermap([9],'BrBg'));
%subplot(2,2,3)
hold on
%plot(nanmean(VPD_BINS(i,:)),[(nanmean(GS_BINS(i,:))-nanstd(GS_BINS(i,:)))'].*1000,'r:');
hold on
%plot(nanmean(VPD_BINS(i,:)),[(nanmean(GS_BINS(i,:))+nanstd(GS_BINS(i,:)))'].*1000,'r:');

plot(nanmean(VPD_BINS(i,:)),GS_BINS(i,:).*1000,'r:','Color',map_CB(2,:));
plot(nanmean(VPD_BINS(i1,:)),GS_BINS(i1,:).*1000,'r:','Color',map_CB(8,:));
l1=plot(nanmean(VPD_BINS(i,:)),nanmean(GS_BINS(i,:)).*1000,'-','Color',map_CB(2,:),'LineWidth',3);
hold on
l2=plot(nanmean(VPD_BINS(i1,:)),nanmean(GS_BINS(i1,:)).*1000,'-','Color',map_CB(8,:),'LineWidth',3);

axis square
legend([l1 l2],{'Peatlands','Forests'});
xlim([-0.1 3.5]);

ylabel('g_s [mm s^{-1}]');
xlabel('VPD [kPa]');
set(gca,'FontSize',16);

%% calculate if gs per VPD is significantly different
% randomly compare peatland and forest gs for bins, then take probability
% from binary distribution
for p=1:40;
    x=GS_BINS(i,p);
    y=GS_BINS(i1,p);
    %p
    x(~isnan(x));
    y(~isnan(y));
    x1=x(~isnan(x));
    y1=y(~isnan(y));
    %pause()
    clc
    for k=1:1000;
       index1 = randsample(1:1:length(x1),1,'true');
       index2 = randsample(1:1:length(y1),1,'true');
       % test if peatland > forest
       PR_TEST_GS(k)=x1(index1)>y1(index2);
       
       %index1 = randsample(1:1:length(x1),length(x1),'true');
       %index2 = randsample(1:1:length(y1),length(y1),'true');
       % test if peatland > forest
       %PR_TEST(k)=mean(x1(index1))>mean(y1(index2));
    end
    PROB_DIFF_GS(p)=sum(PR_TEST_GS)./sum(~isnan(PR_TEST_GS));
%     hist(PR_TEST_GS)
%     pause(.2)
    clear PR_TEST_GS
    
    if length(x1)>=3
        [pval_VPD_GS(p), observeddifference, effectsize] = permutationTest(x1,y1,99);
    elseif length(x1)<3
        pval_VPD_GS(p)=NaN;
    else
        [pval_VPD_GS(p), observeddifference, effectsize] = permutationTest(x1,y1,99,'exact',1);
    end
    clear x y x1 y1
end
%%
r = randi([1 length(i)],1,1000);
r1 = randi([1 length(i1)],1,1000);
for k = 1:1000;
    DeltaGS(:,k)=GS_BINS(i(r(k)),:)-GS_BINS(i1(r1(k)),:);
end
subplot(2,2,4)
%plotshaded(nanmean(VPD_BINS(i,:)),prctile(DeltaET',[15.87 84.14]),'r');
%plot(nanmean(VPD_BINS(i,:)),prctile(DeltaGS',[25]).*1000,'r:');
hold on
%plot(nanmean(VPD_BINS(i,:)),prctile(DeltaGS',[75]).*1000,'r:');
l1=plot(nanmean(VPD_BINS(i,:)),nanmean(DeltaGS').*1000,'-','LineWidth',3);
DET=nanmean(DeltaGS');
%l2=plot(nanmean(VPD_BINS(i,pval_VPD<=0.05)),DET(pval_VPD<=0.05),'.','MarkerSize',14);
%l3=plot(nanmean(VPD_BINS(i,pval_VPD>0.05)),DET(pval_VPD>0.05),'x','MarkerSize',14);
%legend([l2 l3],{'p < 0.05','p > 0.05'});
scatter(nanmean(VPD_BINS),nanmean(DeltaGS').*1000,35,PROB_DIFF_GS,'filled');
cbr=colorbar;
caxis([0.5 0.9]);
set(cbr,'YTickLabel',(0.5:0.1:0.9).*100)
set(get(cbr,'ylabel'),'String', {'Probability of g_{s {PTL}} > g_{s {FOR}}[%]'},'FontSize',14);

ylabel('\Delta g_s [mm s^{-1}]');
xlabel('VPD [kPa]');
set(gca,'FontSize',16);

hold on

axis square
xlim([-0.1 3.5]);
%% calculate if LE per VPD is significantly different
i=find(strcmp(PFT_NUM,'WET'));
i1=find(~strcmp(PFT_NUM,'WET'));
% randomly compare peatland and forest ET for bins, then take probability
% from binary distribution
for p=1:40;
    x=LE_BINS(i,p);
    y=LE_BINS(i1,p);
    p
    x(~isnan(x));
    y(~isnan(y));
    x1=x(~isnan(x));
    y1=y(~isnan(y));
    %pause()
    clc
    for k=1:1000;
       index1 = randsample(1:1:length(x1),1,'true');
       index2 = randsample(1:1:length(y1),1,'true');
       % test if peatland > forest
       PR_TEST(k)=x1(index1)>y1(index2);
       
       %index1 = randsample(1:1:length(x1),length(x1),'true');
       %index2 = randsample(1:1:length(y1),length(y1),'true');
       % test if peatland > forest
       %PR_TEST(k)=mean(x1(index1))>mean(y1(index2));
    end
    PROB_DIFF(p)=sum(PR_TEST)./sum(~isnan(PR_TEST));
    clear PR_TEST
    
    if length(x1)>=3
        [pval_VPD(p), observeddifference, effectsize] = permutationTest(x1,y1,99);
    elseif length(x1)<3
        pval_VPD(p)=NaN;
    else
        [pval_VPD(p), observeddifference, effectsize] = permutationTest(x1,y1,99,'exact',1);
    end
    clear x y x1 y1
end

%%
%figure,
subplot(1,2,1);
% for k = 1:length(VPD_BINS(:,1));
%     if strcmp(PFT_NUM{k},'WET');
%         plot(VPD_BINS(k,:)',LE_BINS(k,:)','-','Color',[0.5 0.5 0.5],'LineWidth',2)
%         hold on
%     else
%        plot(VPD_BINS(k,:)',LE_BINS(k,:)','-','Color',[0.8 0.8 0.8],'LineWidth',2)
%        hold on
%     end
% end

%plot(VPD_BINS_SIT(i1,:)',LE_BINS_SIT(i1,:)','-','Color',[0.65 0.65 0.65],'LineWidth',2);
hold on
% for k = 1:length(i1);
%     if ~isnan(PRE_NUM(i1(k)))
%         plot(VPD_BINS_SIT(i1(k),:)',LE_BINS_SIT(i1(k),:)','-','Color',[PRE_NUM(i1(k))/max(PRE_NUM)-0.1 PRE_NUM(i1(k))/max(PRE_NUM)-0.1 PRE_NUM(i1(k))/max(PRE_NUM)-0.1],'LineWidth',2);
%         hold on
%     end
% end

%plot(VPD_BINS_SIT(i,:)',LE_BINS_SIT(i,:)','-','Color',[0.85 0.85 0.85],'LineWidth',2);
% for k = 1:length(i);
%     if ~isnan(PRE_NUM(i(k)))
%         plot(VPD_BINS_SIT(i(k),:)',LE_BINS_SIT(i(k),:)','-','Color',[PRE_NUM(i(k))/max(PRE_NUM)-0.1 0.5 0.5],'LineWidth',2);
%         hold on
%     end
% end
%plotshaded(nanmean(VPD_BINS(i,:)),[(nanmean(LE_BINS(i,:))-nanstd(LE_BINS(i,:)))' (nanmean(LE_BINS(i,:))+nanstd(LE_BINS(i,:)))'],'r');
plot(nanmean(VPD_BINS(i,:)),[(nanmean(LE_BINS(i,:))-nanstd(LE_BINS(i,:)))'],':','Color',map_CB(2,:));
% (nanmean(LE_BINS(i,:))+nanstd(LE_BINS(i,:)))'
hold on
plot(nanmean(VPD_BINS(i,:)),[(nanmean(LE_BINS(i,:))+nanstd(LE_BINS(i,:)))'],':','Color',map_CB(2,:));
%plotshaded(nanmean(VPD_BINS(find(~strcmp(PFT_NUM,'WET')),:)),[(nanmean(LE_BINS(find(~strcmp(PFT_NUM,'WET')),:))-nanstd(LE_BINS(find(~strcmp(PFT_NUM,'WET')),:)))' (nanmean(LE_BINS(find(~strcmp(PFT_NUM,'WET')),:))+nanstd(LE_BINS(find(~strcmp(PFT_NUM,'WET')),:)))'],'r');
%plotshaded(nanmean(VPD_BINS(i1,:)),[(nanmean(LE_BINS(i1,:))-nanstd(LE_BINS(i1,:)))' (nanmean(LE_BINS(i1,:))+nanstd(LE_BINS(i1,:)))'],'r');
plot(nanmean(VPD_BINS(i1,:)),[(nanmean(LE_BINS(i1,:))+nanstd(LE_BINS(i1,:)))'],':','Color',map_CB(8,:));
plot(nanmean(VPD_BINS(i1,:)),[(nanmean(LE_BINS(i1,:))-nanstd(LE_BINS(i1,:)))'],':','Color',map_CB(8,:));
l1=plot(nanmean(VPD_BINS(i,:)),nanmean(LE_BINS(i,:)),'-','Color',map_CB(2,:),'LineWidth',3);
%plot(nanmean(VPD_BINS(find(~strcmp(PFT_NUM,'WET')),:)),nanmean(LE_BINS(find(~strcmp(PFT_NUM,'WET')),:)),'-');
l2=plot(nanmean(VPD_BINS(i1,:)),nanmean(LE_BINS(i1,:)),'-','Color',map_CB(8,:),'LineWidth',3);
axis square
legend([l1 l2],{'Peatlands','Forests'});
xlim([-0.1 3.5]);
%ylim([0 275]);
ylabel('ET [mm hr^{-1}]');
xlabel('VPD [kPa]');
set(gca,'FontSize',16);
subplot(1,2,2);
r = randi([1 length(i)],1,1000);
r1 = randi([1 length(i1)],1,1000);
for k = 1:1000;
    DeltaET(:,k)=LE_BINS(i(r(k)),:)-LE_BINS(i1(r1(k)),:);
end
%plotshaded(nanmean(VPD_BINS(i,:)),prctile(DeltaET',[15.87 84.14]),'r');
plot(nanmean(VPD_BINS(i,:)),prctile(DeltaET',[25]),'r:');
hold on
plot(nanmean(VPD_BINS(i,:)),prctile(DeltaET',[75]),'r:');
l1=plot(nanmean(VPD_BINS(i,:)),nanmean(DeltaET'),'-','LineWidth',3);
DET=nanmean(DeltaET');
%l2=plot(nanmean(VPD_BINS(i,pval_VPD<=0.05)),DET(pval_VPD<=0.05),'.','MarkerSize',14);
%l3=plot(nanmean(VPD_BINS(i,pval_VPD>0.05)),DET(pval_VPD>0.05),'x','MarkerSize',14);
%legend([l2 l3],{'p < 0.05','p > 0.05'});
scatter(nanmean(VPD_BINS),nanmean(DeltaET'),35,PROB_DIFF,'filled');
cbr=colorbar;
caxis([0.5 0.9]);
set(cbr,'YTickLabel',(0.5:0.1:0.9).*100)
set(get(cbr,'ylabel'),'String', {'Probability of ET_{PTL} > ET_{FOR} [%]'},'FontSize',14);

ylabel('\DeltaET [mm m^{-2} hr^{-1}]');
xlabel('VPD [kPa]');
set(gca,'FontSize',16);

hold on

axis square
xlim([-0.1 3.5]);
%ylim([-0.02 0.23]);
%%
load('D:\Mac_PDF_Boreal\LAKE\LAKE_30_KVJ.mat');
load('D:\Mac_PDF_Boreal\LAKE\LAKE_30_MER.mat');
load('D:\Mac_PDF_Boreal\LAKE\LAKE_30_PNP.mat');
load('D:\Mac_PDF_Boreal\LAKE\LAKE_30_ST1.mat');
load('D:\Mac_PDF_Boreal\LAKE\LAKE_30_VAR.mat');
load('D:\Mac_PDF_Boreal\LAKE\LAKE_30_RES.mat');

subplot(1,2,1);
VPD_LAKE=vertcat(VPD_BINS_KVJ,VPD_BINS_MER,VPD_BINS_PNP,VPD_BINS_ST1,VPD_BINS_VAR,VPD_BINS_RES);
LE_LAKE=vertcat(LE_BINS_KVJ,LE_BINS_MER,LE_BINS_PNP,LE_BINS_ST1,LE_BINS_VAR,LE_BINS_RES);
%plotshaded(nanmean(VPD_LAKE),[(nanmean(LE_LAKE)-nanstd(LE_LAKE))' (nanmean(LE_LAKE)+nanstd(LE_LAKE))'],'r');
plot(nanmean(VPD_LAKE),nanmean(LE_LAKE),'r');
hold on
for k=1:6;
    plot(VPD_LAKE(k,:),LE_LAKE(k,:),'r');
end
%% read monthly data
clearvars -except SITE_NAME PFT SITE TC_ALL SITE_NAME_TC
lng=90;
LE_MNT_PTL=[];
LE_MNT_FOR=[];
LE_MNT_PFR=[];

RN_MNT_PTL=[];
RN_MNT_FOR=[];
RN_MNT_PFR=[];

SWIN_MNT_PTL=[];
SWIN_MNT_FOR=[];

P_MNT_PTL=[];
P_MNT_FOR=[];

VPD_MNT_PTL=[];
VPD_MNT_FOR=[];
VPD_MNT_PFR=[];

VPD_MNT_MID_PTL=[];
VPD_MNT_MID_FOR=[];

SITE_PTL=[];
SITE_FOR=[];
ANN_FOR=[];
MNT_FOR=[];
LAT_FOR=[];
LON_FOR=[];

ANN_PTL=[];
MNT_PTL=[];
LAT_PTL=[];
LON_PTL=[];
PFT_PTL=[];
PFT_FOR=[];

ANN_PFR=[];
MNT_PFR=[];


EVI_MNT_PTL=[];
EVI_MNT_FOR=[];

LAI_MNT_PTL=[];
LAI_MNT_FOR=[];

dEVI_MNT_PTL=[];
dEVI_MNT_FOR=[];

dLAI_MNT_PTL=[];
dLAI_MNT_FOR=[];
PFT_U=unique(PFT);
for k = 1:lng;
    PFT_NUM(k)=find(strcmp(PFT{k},PFT_U));
end
lng=length(SITE_NAME);
load('C:\Users\Manu_MAC\Mega\GWF_Postdoc\PAN_BOREAL\DATA\ESA_CCI_SM_SITES.mat');
%load('C:\Users\Manu_MAC\Mega\GWF_Postdoc\PAN_BOREAL\DATA\EVI_MN_SIT.mat');
load('D:\ET_SYNTHESIS\DATA\LAI_BOR_SITES\LAI_MN_BOR.mat');

%%
%{
EVI_MN_MEAN=NaN(90,12);
dEVI_MN=NaN(size(EVI_MN));
for k=1:90;
    for n=1:12;
        EVI_MN_MEAN(k,n)=nanmean(EVI_MN(k,MN_SIT==n));
        dEVI_MN(k,n:12:end)=EVI_MN(k,n:12:end)-EVI_MN_MEAN(k,n);
    end
end
%}
LAI_MN_MEAN=NaN(90,12);
dLAI_MN=NaN(size(LAI_MN));
for k=1:90;
    for n=1:12;
        LAI_MN_MEAN(k,n)=nanmean(LAI_MN(k,MN_SIT==n));
        dLAI_MN(k,n:12:end)=LAI_MN(k,n:12:end)-LAI_MN_MEAN(k,n);
    end
end

%%
load('D:\ET_SYNTHESIS\TC_SITES.mat');

%% get monthly data and location
%cd('/Users/manuelhelbig/Desktop/Mac_PDF_Boreal/ANN/');
s=0;
t=0;
n=0;
for k = 1:lng; %length(SITE_NAME);
    cd('D:\Mac_PDF_Boreal\MNT\');
    SITE_NAME{k}
    %{
    if strcmp(SITE_NAME{k}(4:end),'YLF')
        SIT_SEL=find(strcmp(SITE_NAME_TC,'SKP'));
    else
        SIT_SEL=find(strcmp(SITE_NAME_TC,SITE_NAME{k}(4:end)));
    end
    %}
    if strcmp(SITE_NAME{k}(4:end),'YLF')
        %SIT_SEL=find(strcmp(SITE_NAME_EVI,'SKP'));
        SIT_SEL=find(strcmp(SITE_NAME_LAI,'SKP'));
    else
        %SIT_SEL=find(strcmp(SITE_NAME_EVI,SITE_NAME{k}(4:end)));
        SIT_SEL=find(strcmp(SITE_NAME_LAI,SITE_NAME{k}(4:end)));
    end
    
    data=importdata(['MNT_' SITE_NAME{k} '.csv']);
    if strcmp(SITE_NAME{k},'FI-ALK') | strcmp(SITE_NAME{k},'FI-KNS')...
            | strcmp(SITE_NAME{k},'FI-LET')| strcmp(SITE_NAME{k},'RU-FYO')...
            | strcmp(SITE_NAME{k},'CA-MAN')| strcmp(SITE_NAME{k},'CA-SCC')
        n=n+1
        ANN_PFR=vertcat(ANN_PFR,data(:,1));
        MNT_PFR=vertcat(MNT_PFR,data(:,2));
        LE_MNT_PFR=vertcat(LE_MNT_PFR,data(:,3));
        VPD_MNT_PFR=vertcat(VPD_MNT_PFR,data(:,4));
        RN_MNT_PFR=vertcat(RN_MNT_PFR,data(:,6));
    end
    if nanmean(data(:,9))==1 % TC_SIT(SIT_SEL)<=40; %
        s=s+1
%         subplot(1,2,1)
%         plot(data(:,4),data(:,2),'.','MarkerSize',24);
%         hold on
%         axis square
%         pause()
        ANN_PTL=vertcat(ANN_PTL,data(:,1));
        MNT_PTL=vertcat(MNT_PTL,data(:,2));
        LE_MNT_PTL=vertcat(LE_MNT_PTL,data(:,3));
        for w=1:length(data(:,1));
            sel=find(YR_SIT==data(w,1) & MN_SIT==data(w,2));
            %{
            if isempty(sel)
                EVI_MNT_PTL1(w)=NaN;
                dEVI_MNT_PTL1(w)=NaN;
            else
                EVI_MNT_PTL1(w)=EVI_MN(SIT_SEL,sel);
                dEVI_MNT_PTL1(w)=dEVI_MN(SIT_SEL,sel);
            end
            %}
            
            if isempty(sel)
                LAI_MNT_PTL1(w)=NaN;
                dLAI_MNT_PTL1(w)=NaN;
            else
                LAI_MNT_PTL1(w)=LAI_MN(SIT_SEL,sel);
                dLAI_MNT_PTL1(w)=dLAI_MN(SIT_SEL,sel);
            end
            
            clear sel
        end
%         EVI_MNT_PTL=vertcat(EVI_MNT_PTL,EVI_MNT_PTL1');
%         dEVI_MNT_PTL=vertcat(dEVI_MNT_PTL,dEVI_MNT_PTL1');
%         clear EVI_MNT_PTL1 dEVI_MNT_PTL1

        LAI_MNT_PTL=vertcat(LAI_MNT_PTL,LAI_MNT_PTL1');
        dLAI_MNT_PTL=vertcat(dLAI_MNT_PTL,dLAI_MNT_PTL1');
        clear LAI_MNT_PTL1 dLAI_MNT_PTL1

%         data(:,3)
%         pause()
        SWIN_MNT_PTL=vertcat(SWIN_MNT_PTL,data(:,7));
        RN_MNT_PTL=vertcat(RN_MNT_PTL,data(:,6));

        VPD_MNT_PTL=vertcat(VPD_MNT_PTL,data(:,4));
        VPD_MNT_MID_PTL=vertcat(VPD_MNT_MID_PTL,data(:,5));
        
        SITE_PTL=vertcat(SITE_PTL,ones(size(data(:,3))).*s);
        %PFT_PTL=vertcat(PFT_PTL,ones(size(data(:,3))).*PFT_NUM(strcmp(SITE_NAME{k},SITE)));
        
        if strcmp(SITE_NAME{k},'RU-YLF');
            PFT_PTL=vertcat(PFT_PTL,ones(size(data(:,3)))*PFT_NUM(strcmp('RU-SKP',SITE)));
        else
            PFT_PTL=vertcat(PFT_PTL,ones(size(data(:,3)))*PFT_NUM(strcmp(SITE_NAME{k},SITE)));
        end
        sze=data(:,4);
    else
        t=t+1
%         subplot(1,2,2)
%         plot(data(:,4),data(:,2),'.','MarkerSize',24);
%         hold on
%         axis square
%         pause()
        ANN_FOR=vertcat(ANN_FOR,data(:,1));
        MNT_FOR=vertcat(MNT_FOR,data(:,2));
        LE_MNT_FOR=vertcat(LE_MNT_FOR,data(:,3));
        
        for w=1:length(data(:,1));
            sel=find(YR_SIT==data(w,1) & MN_SIT==data(w,2));
            %{
            if isempty(sel)
                EVI_MNT_FOR1(w)=NaN;
                dEVI_MNT_FOR1(w)=NaN;
            else
                EVI_MNT_FOR1(w)=EVI_MN(SIT_SEL,sel);
                dEVI_MNT_FOR1(w)=dEVI_MN(SIT_SEL,sel);
            end
            %}
            
            if isempty(sel)
                LAI_MNT_FOR1(w)=NaN;
                dLAI_MNT_FOR1(w)=NaN;
            else
                LAI_MNT_FOR1(w)=LAI_MN(SIT_SEL,sel);
                dLAI_MNT_FOR1(w)=dLAI_MN(SIT_SEL,sel);
            end
            
            clear sel
        end
%         EVI_MNT_FOR=vertcat(EVI_MNT_FOR,EVI_MNT_FOR1');
%         dEVI_MNT_FOR=vertcat(dEVI_MNT_FOR,dEVI_MNT_FOR1');
%         clear EVI_MNT_FOR1 dEVI_MNT_FOR1

        LAI_MNT_FOR=vertcat(LAI_MNT_FOR,LAI_MNT_FOR1');
        dLAI_MNT_FOR=vertcat(dLAI_MNT_FOR,dLAI_MNT_FOR1');
        clear LAI_MNT_FOR1 dLAI_MNT_FOR1

        
        % data(:,3)
        %pause()
        SWIN_MNT_FOR=vertcat(SWIN_MNT_FOR,data(:,7));
        RN_MNT_FOR=vertcat(RN_MNT_FOR,data(:,6));

        VPD_MNT_FOR=vertcat(VPD_MNT_FOR,data(:,4));
        VPD_MNT_MID_FOR=vertcat(VPD_MNT_MID_FOR,data(:,5));
        SITE_FOR=vertcat(SITE_FOR,ones(size(data(:,3))).*t);
        if strcmp(SITE_NAME{k},'RU-YLF');
            PFT_FOR=vertcat(PFT_FOR,ones(size(data(:,3)))*PFT_NUM(strcmp('RU-SKP',SITE)));
        else
            PFT_FOR=vertcat(PFT_FOR,ones(size(data(:,3)))*PFT_NUM(strcmp(SITE_NAME{k},SITE)));
        end
        
        sze=data(:,4);
    end
    %pause()
    
    cd('D:\Mac_PDF_Boreal\ANN\');
        SITE_NAME{k};
        data=importdata(['ANN_' SITE_NAME{k} '.csv']);
        if nanmean(data(:,8))==1 % TC_SIT(SIT_SEL)<=40; %
            LAT_PTL=vertcat(LAT_PTL,data(1,10).*ones(size(sze)));
            LON_PTL=vertcat(LON_PTL,data(1,11).*ones(size(sze)));

        else
            LAT_FOR=vertcat(LAT_FOR,data(1,10).*ones(size(sze)));
            LON_FOR=vertcat(LON_FOR,data(1,11).*ones(size(sze)));
        end
    
        clear sze
end

%% get precip and PET data for individual months
%cd('/Users/manuelhelbig/Desktop/CRU/');
cd('D:\CRU\');
path = '/Users/manuelhelbig/Desktop/CRU/';
path = 'D:\CRU\';
finfo = dir([path '*.nc']);
for k = 1:5;
    if k == 1;
    lat = ncread([path,finfo(k).name],'lat');
    lon = ncread([path,finfo(k).name],'lon');
    end
    if k <= 4;
        PRE(:,:,k*120-119:k*120) = ncread([path,finfo(k*4-2).name],'pre');
    else
        PRE(:,:,k*120-119:k*120-119+83) = ncread([path,finfo(k*4-2).name],'pre');
    end
    
    if k <= 4;
        VAP(:,:,k*120-119:k*120) = ncread([path,finfo(k*4).name],'vap')./10;
    else
        VAP(:,:,k*120-119:k*120-119+83) = ncread([path,finfo(k*4).name],'vap')./10;
    end
    
    
    if k <= 4;
        PET(:,:,k*120-119:k*120) = ncread([path,finfo(k*4-3).name],'pet');
    else
        PET(:,:,k*120-119:k*120-119+83) = ncread([path,finfo(k*4-3).name],'pet');
    end

end
%%
path='C:\Users\Manu_MAC\Documents\CRU\TMX\';
finfo = dir([path '*.nc']);

for k=1:5;
    if k <= 4;
        TMX(:,:,k*120-119:k*120) = ncread([path,finfo(k).name],'tmx');
    else
        TMX(:,:,k*120-119:k*120-119+83) = ncread([path,finfo(k).name],'tmx');
    end
end

MM=[];
YR=[];
for k=1:564/12;
    MM=horzcat(MM,1:12);
    YR=horzcat(YR,ones(1,12).*(1970+k));
end
MM_D=[31 28 31 30 31 30 31 31 30 31 30 31];
for k = 1:564;
    PET(:,:,k)=PET(:,:,k).*MM_D(MM(k));
end

a = 0.611;                       %kPa
b = 17.502;                      %unitless
c = 240.97;                      %C
% ------ Constants needed to convert pressure to concentration or density 
Mw=18;Ma=29;                     %Molar weights of water and air
RR=8.314;                        %universal gas constant (J/(mol K))
%Uses the August-Roche-Magnus formula for solving the Clausius-Clapeyron relation
for k=1:564;
    %es_min(:,:,k)= a .* exp((b.*TMN(:,:,k))./(TMN(:,:,k) + c));   %KPa - staturation vapor pressure
    es_max(:,:,k)= a .* exp((b.*TMX(:,:,k))./(TMX(:,:,k) + c));   %KPa - staturation vapor pressure
end
clear TMX
for s = 1:564;
    VPD_MAX(:,:,s)=es_max(:,:,s)-VAP(:,:,s);
end
%%
for k=1:length(LAT_PTL);
    lon_i=LON_PTL(k);
    lat_i=LAT_PTL(k);
    i1=find((lon-lon_i).^2==min((lon-lon_i).^2),1,'first');
    i2=find((lat-lat_i).^2==min((lat-lat_i).^2),1,'first');
    i=find(YR==ANN_PTL(k) & MM == MNT_PTL(k));
    if ~isempty(i)
        P_MNT_PTL(k)=PRE(i1,i2,i);
        PET_MNT_PTL(k)=PET(i1,i2,i);
        VPD_M_MNT_PTL(k)=VPD_MAX(i1,i2,i);
    else
        P_MNT_PTL(k)=NaN;
        PET_MNT_PTL(k)=NaN;
        VPD_M_MNT_PTL(k)=NaN;
    end
end

%%
for k=1:length(LAT_FOR);
    lon_i=LON_FOR(k);
    lat_i=LAT_FOR(k);
    i1=find((lon-lon_i).^2==min((lon-lon_i).^2),1,'first');
    i2=find((lat-lat_i).^2==min((lat-lat_i).^2),1,'first');
    i=find(YR==ANN_FOR(k) & MM == MNT_FOR(k));
    if ~isempty(i)
        P_MNT_FOR(k)=PRE(i1,i2,i);
        PET_MNT_FOR(k)=PET(i1,i2,i);
        VPD_M_MNT_FOR(k)=VPD_MAX(i1,i2,i);
    else
        P_MNT_FOR(k)=NaN;
        PET_MNT_FOR(k)=NaN;
        VPD_M_MNT_FOR(k)=NaN;
    end
end
%%
path='D:\ET_SYNTHESIS\DATA\Mean_Seasonal_LAI_1653\Mean_Seasonal_LAI_1653\data\LAI_mean_monthly_1981-2015.nc4';
LON_LAI=ncread(path,'lon');
LAT_LAI=ncread(path,'lat');
AVH_LAI=ncread(path,'LAI');
TIME=ncread(path,'time_bnds');
for k=1:length(LAT_PTL);
    lon_i=LON_PTL(k);
    lat_i=LAT_PTL(k);
    i1=find((LON_LAI-lon_i).^2==min((LON_LAI-lon_i).^2),1,'first');
    i2=find((LAT_LAI-lat_i).^2==min((LAT_LAI-lat_i).^2),1,'first');
    
    x=AVH_LAI(i1,i2,12);
    x(x==-9999)=NaN;
    LAI_A_MNT_PTL(k)=nanmean(x);
    clear x
end

for k=1:length(LAT_FOR);
    lon_i=LON_FOR(k);
    lat_i=LAT_FOR(k);
    i1=find((LON_LAI-lon_i).^2==min((LON_LAI-lon_i).^2),1,'first');
    i2=find((LAT_LAI-lat_i).^2==min((LAT_LAI-lat_i).^2),1,'first');
    x=AVH_LAI(i1,i2,12);
    x(x==-9999)=NaN;
    LAI_A_MNT_FOR(k)=nanmean(x);
    clear x
end
%% derive monthly P-PET based on hydrologic year
P_PET=NaN(size(PET));
for k=1:564;
    if k<5 | MM(k)>11 | MM(k)<4;
        P_PET(:,:,k)=NaN;
    else
        if MM(k)==5;
           P_PET(:,:,k)=PRE(:,:,k)-PET(:,:,k);
        else
           P_PET(:,:,k)=P_PET(:,:,k-1)+(PRE(:,:,k)-PET(:,:,k)); 
        end
    end
end

%% derive monthly P-PET based on hydrologic year
%{
P_PET=NaN(size(PET));
for k=1:564;
   P_PET(:,:,k)=PRE(:,:,k)-PET(:,:,k);
end
%}
%% get monthly P-PET
for k=1:length(LAT_PTL);
    lon_i=LON_PTL(k);
    lat_i=LAT_PTL(k);
    i1=find((lon-lon_i).^2==min((lon-lon_i).^2),1,'first');
    i2=find((lat-lat_i).^2==min((lat-lat_i).^2),1,'first');
    if ANN_PTL(k)==2018;
        P_PET_PTL(k)=NaN;
    else
        ind=find(YR==ANN_PTL(k) & MM==MNT_PTL(k));
        P_PET_PTL(k)=P_PET(i1,i2,ind);
    end
    P_PET_SIT_PTL(k)=nanmean(PRE(i1,i2,YR>=2000)-PET(i1,i2,YR>=2000)).*12;
end    

for k=1:length(LAT_FOR);
    lon_i=LON_FOR(k);
    lat_i=LAT_FOR(k);
    i1=find((lon-lon_i).^2==min((lon-lon_i).^2),1,'first');
    i2=find((lat-lat_i).^2==min((lat-lat_i).^2),1,'first');
    if ANN_FOR(k)==2018 | ANN_FOR(k)==2019;
        P_PET_FOR(k)=NaN;
    else
        ind=find(YR==ANN_FOR(k) & MM==MNT_FOR(k));
        P_PET_FOR(k)=P_PET(i1,i2,ind);
    end
    P_PET_SIT_FOR(k)=nanmean(PRE(i1,i2,YR>=2000)-PET(i1,i2,YR>=2000)).*12;
end    

%%
for k = 1:564/12;
    PRE_ANN(:,:,k)=nansum(PRE(:,:,YR==1970+k),3);
    PRE_GS(:,:,k)=nansum(PRE(:,:,(YR==1970+k & (MM>=5 & MM<=9))),3);
    PET_ANN(:,:,k)=nansum(PET(:,:,YR==1970+k),3);
    PET_GS(:,:,k)=nansum(PET(:,:,(YR==1970+k & (MM>=5 & MM<=9))),3);
end
YY=1971:2017;
%%
for k=1:5;
    subplot(2,3,k);
    edges=0.05:0.05:0.95;
    %hist(LE_MNT_PTL(MNT_PTL==k+4)./PET_MNT_PTL(MNT_PTL==k+4)')
    l1=histogram(LE_MNT_PTL(MNT_PTL==k+4)./PET_MNT_PTL(MNT_PTL==k+4)','BinEdges',edges, 'Normalization','probability', 'DisplayStyle','stairs')
    PET_FR(1,k)=nanmedian(LE_MNT_PTL(MNT_PTL==k+4)./PET_MNT_PTL(MNT_PTL==k+4)');
    hold on
    l2=histogram(LE_MNT_FOR(MNT_FOR==k+4)./PET_MNT_FOR(MNT_FOR==k+4)','BinEdges',edges, 'Normalization','probability', 'DisplayStyle','stairs')
    %hist(LE_MNT_FOR(MNT_FOR==k+4)./PET_MNT_FOR(MNT_FOR==k+4)')
    PET_FR(2,k)=nanmedian(LE_MNT_FOR(MNT_FOR==k+4)./PET_MNT_FOR(MNT_FOR==k+4)');
    xlim([0 1])
end
%% plot VPD-ET sensitivity for all sites with individual slopes
i1=5;
i2=9;
x_FOR=VPD_MNT_MID_FOR(MNT_FOR>=i1 & MNT_FOR<=i2 & PFT_FOR<=6);
y_FOR=LE_MNT_FOR(MNT_FOR>=i1 & MNT_FOR<=i2 & PFT_FOR<=6);
z_FOR=P_MNT_FOR(MNT_FOR>=i1 & MNT_FOR<=i2 & PFT_FOR<=6);
r_FOR=RN_MNT_FOR(MNT_FOR>=i1 & MNT_FOR<=i2 & PFT_FOR<=6);
v_FOR=SM_FOR(MNT_FOR>=i1 & MNT_FOR<=i2 & PFT_FOR<=6);
w_FOR=LAI_MNT_FOR(MNT_FOR>=i1 & MNT_FOR<=i2 & PFT_FOR<=6);
p_FOR=LAI_A_MNT_FOR(MNT_FOR>=i1 & MNT_FOR<=i2 & PFT_FOR<=6);
%i=find(~isnan(x_FOR) & ~isnan(y_FOR));% & ~isnan(z_FOR));
i=find(~isnan(x_FOR) & ~isnan(y_FOR));
x_FOR=x_FOR(i);
y_FOR=y_FOR(i);
r_FOR=r_FOR(i);
z_FOR=z_FOR(i);
w_FOR=w_FOR(i);
v_FOR=v_FOR(i);
p_FOR=p_FOR(i);
%%
x_PFR=VPD_MNT_MID_PFR(MNT_PFR>=i1 & MNT_PFR<=i2);
y_PFR=LE_MNT_PFR(MNT_PFR>=i1 & MNT_PFR<=i2);
r_PFR=RN_MNT_PFR(MNT_PFR>=i1 & MNT_PFR<=i2);
i=find(~isnan(x_PFR) & ~isnan(y_PFR));
x_PFR=x_PFR(i);
y_PFR=y_PFR(i);
r_PFR=r_PFR(i);

%%
x_PTL=VPD_MNT_MID_PTL(MNT_PTL>=i1 & MNT_PTL<=i2);
y_PTL=LE_MNT_PTL(MNT_PTL>=i1 & MNT_PTL<=i2);
z_PTL=P_MNT_PTL(MNT_PTL>=i1 & MNT_PTL<=i2);
r_PTL=RN_MNT_PTL(MNT_PTL>=i1 & MNT_PTL<=i2);
w_PTL=LAI_MNT_PTL(MNT_PTL>=i1 & MNT_PTL<=i2);
v_PTL=SM_PTL(MNT_PTL>=i1 & MNT_PTL<=i2);
i=find(~isnan(x_PTL) & ~isnan(y_PTL)); % & ~isnan(z_PTL));
%i=find(~isnan(x_PTL) & ~isnan(y_PTL) & ~isnan(w_PTL));
x_PTL=x_PTL(i);
y_PTL=y_PTL(i);
r_PTL=r_PTL(i);
z_PTL=z_PTL(i);
w_PTL=w_PTL(i);
v_PTL=v_PTL(i);
%%
%{
figure,
ind=prctile(unique(z_FOR),1:10:90);

VPD_RESP=@(par,VPD) ((-par(1).*log(VPD)+par(2)).*VPD);
    opts = statset('nlinfit');
    opts.RobustWgtFun = 'bisquare';

clear ind
ind=prctile(unique(z_PTL),1:10:90);
for k=1:10;
    x=x_PTL(z_PTL>=ind(k) & z_PTL<ind(k+1));
    y=y_PTL(z_PTL>=ind(k) & z_PTL<ind(k+1));
    plot(x,y,'.','MarkerSize',12,'Color',[k/10 k/10 k/10]);
    lin_par(:,k)=polyfit(x,y,1);
    par(:,k) = nlinfit(x,y,VPD_RESP,[23 61],opts);
    hold on
    plot([min(x):0.1:max(x)],lin_par(1,k).*[min(x):0.1:max(x)]+lin_par(2,k),'-','LineWidth',2,'Color',[k/10 k/10 k/10])
    plot([min(x):0.1:max(x)],VPD_RESP(par(:,k),[min(x):0.1:max(x)]),'-','LineWidth',2,'Color',[k/10 k/10 k/10])
end

for k=1:10;
    x=x_FOR(z_FOR>=ind(k) & z_FOR<ind(k+1));
    y=y_FOR(z_FOR>=ind(k) & z_FOR<ind(k+1));
    plot(x,y,'.','MarkerSize',12,'Color',[k/10 k/10 k/10]);
    lin_par(:,k)=polyfit(x,y,1);
    par(:,k) = nlinfit(x,y,VPD_RESP,[23 61],opts);
    hold on
    plot([min(x):0.1:max(x)],lin_par(1,k).*[min(x):0.1:max(x)]+lin_par(2,k),'-','LineWidth',2,'Color',[k/10 k/10 k/10])
    plot([min(x):0.1:max(x)],VPD_RESP(par(:,k),[min(x):0.1:max(x)]),'-','LineWidth',2,'Color',[k/10 k/10 k/10])
end
%}
VPD_RESP=@(par,VPD) ((1./(1+(VPD./par(1))).*par(2).*VPD));
% see Oren 1999 (accounts for decrease in stomatal conductance with VPD)
VPD_RESP=@(par,VPD) ((-par(1).*log(VPD)+par(2)).*VPD);
    opts = statset('nlinfit');
    opts.RobustWgtFun = 'bisquare';
%VPD_RESP=@(par,VPD) (((-par(1).*log(VPD(:,1))+(par(2))).*(VPD(:,2)./max(VPD(:,2)))).*VPD(:,1));
    
figure,
subplot(1,2,1)
thr=1000;
for k=1;
    if k==1;
        x=x_PTL;
        y=y_PTL;
    else
        x=x_PFR;
        y=y_PFR;
    end
    clear par
    w=w_PTL;
    for s=1:23;
        y_bin_md(s)=nanmedian(y_PTL(x_PTL>(s./10)-0.1 & x_PTL<=(s./10)));
        y_bin_mn(s)=nanmean(y_PTL(x_PTL>(s./10)-0.1 & x_PTL<=(s./10)));
    end
    plot(0.05:0.1:2.25,y_bin_md);
    hold on 
    plot(0.05:0.1:2.25,y_bin_mn);
    %{
    if k == 1;
        x=x_PTL(z_PTL<=thr);
        y=y_PTL(z_PTL<=thr);
        w=w_PTL(z_PTL<=thr);
    else
        x=x_PTL(z_PTL>thr);
        y=y_PTL(z_PTL>thr);
        w=w_PTL(z_PTL>thr);
    end
    %}
    %pr=polyfit(w,y,1);
    %y=y-(pr(1)*w+pr(2));
    clear par unc_sl
    %x=horzcat(x,w);
    for n = 1:1000;

        if n==1;
            %par(:,n) = polyfit(x,y,1);
            %unc_sl(:,n)=par(1).*(min(x):0.1:max(x))+par(2);
            par(:,n) = nlinfit(x,y,VPD_RESP,[23 61],opts);
            %unc_sl(:,n)=VPD_RESP(par(:,n),[min(x):0.1:max(x)]);
        else
            index = randsample(1:1:length(x(:,1)),length(x(:,1)),'true');
            par(:,n) = nlinfit(x(index,:),y(index,:),VPD_RESP,[23 61],opts);
            %unc_sl(:,n)=VPD_RESP(par(:,n),[min(x):0.1:max(x)]);
%         else
%             index = randsample(1:1:length(x),length(x),'true');
%             par(:,n) = polyfit(x(index),y(index),1);
%             unc_sl(:,n)=par(1,n).*(min(x):0.1:max(x))+par(2,n);
        end
    end
    par_PTL{k}=par;
    m_PTL(k)=par(1,1);
    low_CI_PTL(k)=prctile(par(1,:),2.5);
    up_CI_PTL(k)=prctile(par(1,:),97.5);
    b_PTL(k)=par(2,1);
    low_b_CI_PTL(k)=prctile(par(2,:),2.5);
    up_b_CI_PTL(k)=prctile(par(2,:),97.5);
    % pr=polyfit(x,y,1);
    % pr
    l{k}=plot(x,y,'.','MarkerSize',20);
    %l{k}=scatter(x,y,20,w,'filled');
    %caxis([0 0.55]);
    hold on
    PRCT=prctile(x,[5:5:95]);
    PRCT_FIX=0.1:0.05:2;
%     for k=1:20;
%         if k==1;
%             i=find(x<PRCT(k));
%         elseif k==20;
%             i=find(x>=PRCT(k-1));
%         else
%            i=find(x>=PRCT(k-1) & x<PRCT(k)); 
%         end    
%         %i=find(x>k./10 & x<=k./10+0.1);
%         x_PTL_MN(k)=nanmean(x(i));
%         y_PTL_MN(k)=nanmean(y(i));
% 
%     end

%     for k=1:39;
%         i=find(x>=PRCT_FIX(k) & x<PRCT_FIX(k)+0.25);
% 
%         %i=find(x>k./10 & x<=k./10+0.1);
%         x_PTL_MN(k)=nanmean(x(i));
%         y_PTL_MN(k)=nanmean(y(i));
% 
%     end
    %plot(x_PTL_MN,y_PTL_MN,'-','LineWidth',3);
    [R p]=corrcoef(x,y);
    R2_PTL(k,1)=R(1,2)^2;
    lin_par=polyfit(x,y,1);
    %[Err, lin_par] = fit_2D_data(x, y,'no');
    RMSE(k,1)=sqrt(nanmean((y-(lin_par(1).*x+lin_par(2))).^2));
    [R p]=corrcoef(VPD_RESP(par(:,1),x),y);
    hold on
    plot(0:0.1:2,VPD_RESP(par(:,1),(0:0.1:2)),'-');
    %hold on
    %plot(0:125,0:125)
    %xlim([0 125]);
    %ylim([0 125]);
    R2_PTL(k,2)=R(1,2)^2;
    RMSE(k,2)=sqrt(nanmean((y-VPD_RESP(par(:,1),x)).^2));
    pval_PTL=p(1,2);
    ksr_out=ksr(x,y,0.2,length([min(x):0.01:max(x)]));
    %plot(ksr_out.x,ksr_out.f);
    %plotshaded([min(x):0.1:max(x)],prctile(unc_sl',[2.5 97.5]),'r');
    %plot([min(x):0.1:max(x)],lin_par(1).*[min(x):0.1:max(x)]+lin_par(2),'-','LineWidth',2)
   % plot([min(x):0.1:max(x)],VPD_RESP(par(:,1),[min(x):0.1:max(x)]),'-','LineWidth',2)
    axis square
    res_PTL=y_PTL-VPD_RESP(par(:,1),x);
end

%%
x_PTL=VPD_MNT_MID_PTL(MNT_PTL>=i1 & MNT_PTL<=i2);
y_PTL=LE_MNT_PTL(MNT_PTL>=i1 & MNT_PTL<=i2);
z_PTL=SITE_PTL(MNT_PTL>=i1 & MNT_PTL<=i2);
m_PTL= MNT_PTL(MNT_PTL>=i1 & MNT_PTL<=i2);
%w_PTL=P_PET_PTL(MNT_PTL>=i1 & MNT_PTL<=i2)';
%r_PTL=P_MNT_PTL(MNT_PTL>=i1 & MNT_PTL<=i2)';
%w_PTL=(y_PTL./PET_MNT_PTL(MNT_PTL>=i1 & MNT_PTL<=i2)');
%w_PTL=LAI_MNT_PTL(MNT_PTL>=i1 & MNT_PTL<=i2);
t_PTL=SM_PTL(MNT_PTL>=i1 & MNT_PTL<=i2);
i=find(~isnan(x_PTL) & ~isnan(y_PTL));
x_PTL=x_PTL(i);
y_PTL=y_PTL(i);
z_PTL=z_PTL(i);
%r_PTL=r_PTL(i);
%w_PTL=w_PTL(i);
t_PTL=t_PTL(i);
r=0;
%{
for k=1:length(unique(SITE_PTL));
    
    x=x_PTL(z_PTL==k);
    y=y_PTL(z_PTL==k);
    w=w_PTL(z_PTL==k);
    if length(x)>4;
        r=r+1;
        DI_PLT(r)=nanmean(w);
        pr(r,:) = nlinfit(x,y,VPD_RESP,[0.7 130],opts);
        if DI_PLT(r)<0
            plot([min(x):0.1:max(x)],VPD_RESP(pr(r,:),[min(x):0.1:max(x)]),'-','Color',[192/255 192/255 192/255]);
        else
            plot([min(x):0.1:max(x)],VPD_RESP(pr(r,:),[min(x):0.1:max(x)]),'-','Color',[255/255 204/255 204/255]);
        end
    end
end
%}
%PARAM_PLT_SITES=pr;

ylabel('ET_{m} [mm month^{-1}]');
xlabel('VPD_{m} [kPa]');
set(gca,'FontSize',14);
%% get lake data

load('D:\Mac_PDF_Boreal\LAKE\LAKE_MNT_KVJ_ALL.mat');
%load('D:\Mac_PDF_Boreal\LAKE\LAKE_MNT_MER.mat');
load('D:\Mac_PDF_Boreal\LAKE\LAKE_MNT_PNP_ALL.mat');
load('D:\Mac_PDF_Boreal\LAKE\LAKE_MNT_ST1_ALL.mat');
load('D:\Mac_PDF_Boreal\LAKE\LAKE_MNT_VAR_ALL.mat');
load('D:\Mac_PDF_Boreal\LAKE\LAKE_MNT_RES_ALL.mat');
load('D:\Mac_PDF_Boreal\LAKE\LAKE_MNT_TAM_ALL.mat');
load('D:\Mac_PDF_Boreal\LAKE\LAKE_MNT_ERK_ALL.mat');
load('D:\Mac_PDF_Boreal\LAKE\LAKE_MNT_GSL.mat');
load('D:\Mac_PDF_Boreal\LAKE\LAKE_MNT_KOT_ALL.mat');

%Rn_LAKE=horzcat(Rn_MNT_MID_KVJ,Rn_MNT_MID_PNP,Rn_MNT_MID_ST1);
VPD_LAKE=horzcat(VPD_MNT_MID_GSL(MNT_GSL>=i1 & MNT_GSL<=i2),VPD_MNT_MID_KVJ(MNT_KVJ>=i1 & MNT_KVJ<=i2),VPD_MNT_MID_PNP(MNT_PNP>=i1 & MNT_PNP<=i2),VPD_MNT_MID_ST1(MNT_ST1>=i1 & MNT_ST1<=i2),VPD_MNT_MID_VAR(MNT_VAR>=i1 & MNT_VAR<=i2),VPD_MNT_MID_RES(MNT_RES>=i1 & MNT_RES<=i2),VPD_MNT_MID_ERK(MNT_ERK>=i1 & MNT_ERK<=i2),VPD_MNT_MID_TAM(MNT_TAM>=i1 & MNT_TAM<=i2),VPD_MNT_MID_KOT(MNT_KOT>=i1 & MNT_KOT<=i2));
LE_LAKE=horzcat(VPD_MNT_MID_GSL(MNT_GSL>=i1 & MNT_GSL<=i2),LE_MNT_KVJ(MNT_KVJ>=i1 & MNT_KVJ<=i2),LE_MNT_PNP(MNT_PNP>=i1 & MNT_PNP<=i2),LE_MNT_ST1(MNT_ST1>=i1 & MNT_ST1<=i2),LE_MNT_VAR(MNT_VAR>=i1 & MNT_VAR<=i2),LE_MNT_RES(MNT_RES>=i1 & MNT_RES<=i2),LE_MNT_ERK(MNT_ERK>=i1 & MNT_ERK<=i2),LE_MNT_TAM(MNT_TAM>=i1 & MNT_TAM<=i2),LE_MNT_KOT(MNT_KOT>=i1 & MNT_KOT<=i2));
MNT_LAKE=horzcat(VPD_MNT_MID_GSL(MNT_GSL>=i1 & MNT_GSL<=i2),MNT_KVJ(MNT_KVJ>=i1 & MNT_KVJ<=i2),MNT_PNP(MNT_PNP>=i1 & MNT_PNP<=i2),MNT_ST1(MNT_ST1>=i1 & MNT_ST1<=i2),MNT_VAR(MNT_VAR>=i1 & MNT_VAR<=i2),MNT_RES(MNT_RES>=i1 & MNT_RES<=i2),MNT_ERK(MNT_ERK>=i1 & MNT_ERK<=i2),MNT_TAM(MNT_TAM>=i1 & MNT_TAM<=i2),MNT_KOT(MNT_KOT>=i1 & MNT_KOT<=i2));
SITE_LAKE=horzcat(ones(1,length(VPD_MNT_MID_GSL(MNT_GSL>=i1 & MNT_GSL<=i2))),...
    ones(1,length(MNT_KVJ(MNT_KVJ>=i1 & MNT_KVJ<=i2))).*2,...
    ones(1,length(MNT_PNP(MNT_PNP>=i1 & MNT_PNP<=i2))).*3,...
    ones(1,length(MNT_ST1(MNT_ST1>=i1 & MNT_ST1<=i2))).*4,...
    ones(1,length(MNT_VAR(MNT_VAR>=i1 & MNT_VAR<=i2))).*5,...
    ones(1,length(MNT_RES(MNT_RES>=i1 & MNT_RES<=i2))).*6,...
    ones(1,length(MNT_ERK(MNT_ERK>=i1 & MNT_ERK<=i2))).*7,...
    ones(1,length(MNT_TAM(MNT_TAM>=i1 & MNT_TAM<=i2))).*8,...
    ones(1,length(MNT_KOT(MNT_KOT>=i1 & MNT_KOT<=i2))).*9);

x=VPD_LAKE;
y=LE_LAKE;
i=find(~isnan(x) & ~isnan(y));
x=x(i);
y=y(i);
%l{3}=plot(x,y,'.','MarkerSize',20);
plot(x,y,'.','MarkerSize',20);
clear par unc_sl
for n = 1:1000;
        
    if n==1;
        
         
        par(:,n) = polyfit(x,y,1);
        unc_sl(:,n)=par(1).*(min(x):0.1:max(x))+par(2);
    else
        index = randsample(1:1:length(x),length(x),'true');
        par(:,n) = polyfit(x(index),y(index),1);
        unc_sl(:,n)=par(1,n).*(min(x):0.1:max(x))+par(2,n);
    end
end
par_LAKE=par;
low_CI_LAKE(k)=prctile(par(1,:),2.5);
up_CI_LAKE(k)=prctile(par(1,:),97.5);
low_CI_I_LAKE(k)=prctile(par(2,:),2.5);
up_CI_I_LAKE(k)=prctile(par(2,:),97.5);

[R p]=corrcoef(x,y);
R2_LAK=R(1,2)^2;
pval_LAK=p(1,2);
plotshaded([min(x):0.1:max(x)],prctile(unc_sl',[2.5 97.5]),'r');
hold on
plot([min(x):0.1:max(x)],par(1,1).*[min(x):0.1:max(x)]+par(2,1),'-','LineWidth',2)
%legend([l{1} l{2} l{3}],{'PTL: P<PET','PTL: P>PET','LAKE'});
%ylabel('ET_{mnth} [mm]');
%xlabel('VPD_{max, mnth} [kPa]');
%set(gca,'FontSize',14);
%%
subplot(1,2,1)
clear x y w l
for k=1;
    x=x_FOR;
    y=y_FOR;
   % w=w_FOR;
    for s=1:23;
        y_bin_md(s)=nanmedian(y_FOR(x_FOR>(s./10)-0.1 & x_FOR<=(s./10)));
        y_bin_mn(s)=nanmean(y_FOR(x_FOR>(s./10)-0.1 & x_FOR<=(s./10)));
    end
    plot(0.05:0.1:2.25,y_bin_md);
    hold on 
    plot(0.05:0.1:2.25,y_bin_mn);
    %{
    if k==1;
        x=x_FOR(z_FOR<=thr);
        y=y_FOR(z_FOR<=thr);
        w=w_FOR(z_FOR<=thr);
    else
        x=x_FOR(z_FOR>thr);
        y=y_FOR(z_FOR>thr);
        w=w_FOR(z_FOR>thr);
    end
    %}
    %pr=polyfit(w,y,1);
    %y=y-(pr(1)*w+pr(2));
    %VPD_RESP=@(par,VPD) ((1./(1+(VPD./par(1))).*par(2).*VPD));
    %VPD_RESP=@(par,VPD) (par(3).*(VPD(:,2)./0.6).*(-par(1).*log(VPD(:,1))+par(2)).*VPD(:,1));
    %VPD_RESP=@(par,VPD) (((-par(1).*log(VPD(:,1))+(par(2))).*(VPD(:,2)./max(VPD(:,2)))).*VPD(:,1));
    %x=horzcat(x,w);
    %VPD_RESP=@(par,VPD) ((-par(1).*log(VPD)+par(2)).*VPD);
    opts = statset('nlinfit');
    opts.RobustWgtFun = 'bisquare';
    clear par unc_sl
    for n = 1:1000;
        
        if n==1;
            %par(:,n) = polyfit(x,y,1);
            %RMSE(k*2-1)=sqrt(nanmean((y-(par(1,1).*x+par(2,1))).^2));
            %plot(x,(par(1,1).*x+par(2,1)));
            par(:,n) = nlinfit(x,y,VPD_RESP,[70 125],opts);
            %unc_sl(:,n)=VPD_RESP(par(:,n),[min(x):0.1:max(x)]);
            %RMSE(k*2)=sqrt(nanmean((y-VPD_RESP(par(:,n),x)).^2));
        else
            index = randsample(1:1:length(x),length(x),'true');
            par(:,n) = nlinfit(x(index,:),y(index,:),VPD_RESP,[70 125],opts);
            %unc_sl(:,n)=VPD_RESP(par(:,n),[min(x):0.1:max(x)]);
        end
    end
    par_FOR{k}=par;
    m_FOR(k)=par(1,1);
    low_CI_FOR(k)=prctile(par(1,:),2.5);
    up_CI_FOR(k)=prctile(par(1,:),97.5);
    b_FOR(k)=par(2,1);
    low_b_CI_FOR(k)=prctile(par(2,:),2.5);
    up_b_CI_FOR(k)=prctile(par(2,:),97.5);
%     c_FOR(k)=par(3,1);
%     low_c_CI_FOR(k)=prctile(par(3,:),2.5);
%     up_c_CI_FOR(k)=prctile(par(3,:),97.5);
    %pr = nlinfit(x,y,VPD_RESP,[0.7 130],opts);
    %y=y-(pr(1)*w+pr(2));
    %pr=polyfit(log(x),y,1);
    %pr=polyfit(x,y,1);
    %pr
    l{k}=plot(x,y,'.','MarkerSize',20);
    hold on
    %l{k}=scatter(x,y,20,w,'filled');
    %caxis([0 0.55]);
    hold on
    % for k=1:20;
    %     if k==1;
    %         i=find(x<PRCT(k));
    %     elseif k==20;
    %         i=find(x>=PRCT(k-1));
    %     else
    %        i=find(x>=PRCT(k-1) & x<PRCT(k)); 
    %     end    
    %     %i=find(x_FOR>k./10 & x_FOR<=k./10+0.1);
    %     x_FOR_MN(k)=nanmean(x(i));
    %     y_FOR_MN(k)=nanmean(y(i));
    % end
    % 
    % for k=1:39;
    %     i=find(x>=PRCT_FIX(k) & x<PRCT_FIX(k)+0.25);
    %   
    %     %i=find(x>k./10 & x<=k./10+0.1);
    %     x_FOR_MN(k)=nanmean(x(i));
    %     y_FOR_MN(k)=nanmean(y(i));
    %     
    % end
    %plot(x_FOR_MN,y_FOR_MN,'-','LineWidth',3);
    [R p]=corrcoef(x(:,1),y);
    R2_FOR(k,1)=R(1,2)^2;
    lin_par=polyfit(x(:,1),y,1);
    %[Err, lin_par] = fit_2D_data(x, y,'no');
    RMSE_FOR(k,1)=sqrt(nanmean((y-(lin_par(1).*x(:,1)+lin_par(2))).^2));
    [R p]=corrcoef(VPD_RESP(par(:,1),x),y);
%     plot(VPD_RESP(par(:,1),x),y,'o');
%     hold on
%     plot(0:125,0:125)
    axis square
%     xlim([0 125]);
%     ylim([0 125]);
    R2_FOR(k,2)=R(1,2)^2;
    % including EVI results in R2 of 0.28
    RMSE_FOR(k,2)=sqrt(nanmean((y-VPD_RESP(par(:,1),x)).^2));
    pval_FOR=p(1,2);
    %ksr_out=ksr(x,y,0.2,length([min(x):0.01:max(x)]));
    %plot(ksr_out.x,ksr_out.f);
    %plot([min(x):0.1:max(x)],pr(1).*log([min(x):0.1:max(x)])+pr(2),'-','LineWidth',2)
    %plot([min(x):0.1:max(x)],lin_par(1).*[min(x):0.1:max(x)]+lin_par(2),'-','LineWidth',2)
    %plotshaded([min(x):0.1:max(x)],prctile(unc_sl',[2.5 97.5]),'r');
    plot([min(x):0.1:max(x)],VPD_RESP(par(:,1),[min(x):0.1:max(x)]),'-','LineWidth',2)
    %slope_FOR=pr(1);
    %axis square
    res_FOR=y_FOR-VPD_RESP(par(:,1),x);
end
%legend([l{1} l{2}],{'FOR: P<PET','FOR: P>PET'});
set(gca,'FontSize',14);
%%
x_FOR=VPD_MNT_MID_FOR(MNT_FOR>=i1 & MNT_FOR<=i2 & PFT_FOR<=6);
y_FOR=LE_MNT_FOR(MNT_FOR>=i1 & MNT_FOR<=i2 & PFT_FOR<=6);
z_FOR=SITE_FOR(MNT_FOR>=i1 & MNT_FOR<=i2 & PFT_FOR<=6);
m_FOR=MNT_FOR(MNT_FOR>=i1 & MNT_FOR<=i2 & PFT_FOR<=6);
%w_FOR=P_PET_SIT_FOR(MNT_FOR>=i1 & MNT_FOR<=i2)';
r_FOR=P_MNT_FOR(MNT_FOR>=i1 & MNT_FOR<=i2)';
t_FOR=SM_FOR(MNT_FOR>=i1 & MNT_FOR<=i2);
%w_FOR=(y_FOR./PET_MNT_FOR(MNT_FOR>=i1 & MNT_FOR<=i2 & PFT_FOR<=6)');
w_FOR=LAI_MNT_FOR(MNT_FOR>=i1 & MNT_FOR<=i2 & PFT_FOR<=6);
t_FOR=SM_FOR(MNT_FOR>=i1 & MNT_FOR<=i2)';
i=find(~isnan(x_FOR) & ~isnan(y_FOR));
x_FOR=x_FOR(i);
y_FOR=y_FOR(i);
z_FOR=z_FOR(i);
r_FOR=r_FOR(i);
m_FOR=m_FOR(i);
w_FOR=w_FOR(i);
t_FOR=t_FOR(i);
r=0;
%{
for k=1:length(unique(SITE_FOR));
    
    x=x_FOR(z_FOR==k);
    y=y_FOR(z_FOR==k);
    w=w_FOR(z_FOR==k);
    if length(x)>4;
        r=r+1;
        DI_FOR(r)=nanmean(w);
        pr(r,:) = nlinfit(x,y,VPD_RESP,[0.7 130],opts);
        if DI_FOR(r)<0
            plot([min(x):0.1:max(x)],VPD_RESP(pr(r,:),[min(x):0.1:max(x)]),'-','Color',[192/255 192/255 192/255]);
        else
            plot([min(x):0.1:max(x)],VPD_RESP(pr(r,:),[min(x):0.1:max(x)]),'-','Color',[255/255 204/255 204/255]);
        end
    end
end
%}
%PARAM_FOR_SITES=pr;
%%
map=flipud(brewermap([7],'BrBg'));
figure,
subplot(1,2,1)
plot(w_FOR,res_FOR,'.','MarkerSize',16,'Color',map(2,:));
[p,h] = ranksum(w_FOR,res_FOR)
[R pval]=corrcoef(w_FOR(~isnan(w_FOR)&~isnan(res_FOR)),res_FOR(~isnan(w_FOR)&~isnan(res_FOR)));
x=w_FOR(~isnan(w_FOR)&~isnan(res_FOR));
y=res_FOR(~isnan(w_FOR)&~isnan(res_FOR));
par_res_FOR=NaN(2,1000);
unc_sl=NaN(length([min(x):0.05:max(x)]),1000);
for n = 1:1000;

        if n==1;
            par_res_FOR(:,n)=polyfit(x,y,1);
            unc_sl(:,n)=par_res_FOR(1,n).*([min(x):0.05:max(x)])+par_res_FOR(2,n);
        else
            index = randsample(1:1:length(x),length(x),'true');
            par_res_FOR(:,n) = polyfit(x(index),y(index),1);
            unc_sl(:,n)=par_res_FOR(1,n).*([min(x):0.05:max(x)])+par_res_FOR(2,n);
        end
end

R(1,2)^2

hold on
plotshaded([min(x):0.05:max(x)],prctile(unc_sl',[2.5 97.5]),'r');
plot([min(x)-0.05:0.05:max(x)+0.05],par_res_FOR(1).*([min(x)-0.05:0.05:max(x)+0.05])+par_res_FOR(2),'LineWidth',2,'Color',map(1,:));
clear unc_sl par_res_FOR

refline(0,0);
axis square
subplot(1,2,2)
plot(w_PTL,res_PTL,'.','MarkerSize',16,'Color',map(6,:));
hold on
[p,h] = ranksum(w_PTL,res_PTL)
R=corrcoef(w_PTL(~isnan(w_PTL)&~isnan(res_PTL)),res_PTL(~isnan(w_PTL)&~isnan(res_PTL)));

x=w_PTL(~isnan(w_PTL)&~isnan(res_PTL));
y=res_PTL(~isnan(w_PTL)&~isnan(res_PTL));
par_res_PTL=NaN(2,1000);
unc_sl=NaN(length([min(x):0.05:max(x)]),1000);
for n = 1:1000;

        if n==1;
            par_res_PTL(:,n)=polyfit(x,y,1);
            unc_sl(:,n)=par_res_PTL(1,n).*([min(x):0.05:max(x)])+par_res_PTL(2,n);
        else
            index = randsample(1:1:length(x),length(x),'true');
            par_res_PTL(:,n) = polyfit(x(index),y(index),1);
            unc_sl(:,n)=par_res_PTL(1,n).*([min(x):0.05:max(x)])+par_res_PTL(2,n);
        end
end

R(1,2)^2
hold on
plotshaded([min(x):0.05:max(x)],prctile(unc_sl',[2.5 97.5]),'r');
plot([min(x)-0.05:0.05:max(x)+0.05],par_res_PTL(1).*([min(x)-0.05:0.05:max(x)+0.05])+par_res_PTL(2),'LineWidth',2,'Color',map(1,:));
clear unc_sl par_res_PTL
refline(0,0);

%%
figure,
subplot(1,2,1);
scatter(x_FOR,y_FOR,20,w_FOR,'filled')
caxis([0.4 1])
w_FOR(w_FOR>1)=1;
map=flipud(brewermap([3],'BrBg'));
for k=1:3;
    if k==1;
        ind=find(w_FOR<=0.5);
    else
        ind=find(w_FOR<=0.5+0.25*k-0.25 & w_FOR>0.5+0.25*k-0.5);
    end
    hold on
    n_FOR(k)=length(ind);
    for n = 1:1000;
        
        if n==1;
            %par_DI_B(:,n)= nlinfit(x_FOR(ind),y_FOR(ind),VPD_RESP,[0.7 130],opts);
             par_DI_B(:,n)= polyfit(x_FOR(ind),y_FOR(ind),1);
            %unc_sl(:,n)=VPD_RESP(par_DI_B(:,n),[min(x_FOR):0.1:max(x_FOR)]);
            unc_sl(:,n)=par_DI_B(1,n).*([min(x_FOR):0.1:max(x_FOR)])+par_DI_B(2,n);
        else
            x=x_FOR(ind);
            y=y_FOR(ind);
            index = randsample(1:1:length(x),length(x),'true');
            %par_DI_B(:,n)= nlinfit(x(index),y(index),VPD_RESP,[0.7 130],opts);
            par_DI_B(:,n)= polyfit(x(index),y(index),1);
            %unc_sl(:,n)=VPD_RESP(par_DI_B(:,n),[min(x_FOR):0.1:max(x_FOR)]);
            unc_sl(:,n)=par_DI_B(1,n).*([min(x_FOR):0.1:max(x_FOR)])+par_DI_B(2,n);
        end
    end
    
    par_DI(:,k)= par_DI_B(:,1);
    par_DI_CI{k}=prctile(par_DI_B',[2.5 97.5]);
    % par_DI(:,k)= polyfit(x_FOR(ind),y_FOR(ind),1);
    
    plotshaded([min(x_FOR):0.1:max(x_FOR)],prctile(unc_sl',[2.5 97.5]),'r');
    hold on
     plot([min(x_FOR):0.1:max(x_FOR)],par_DI(1,k).*([min(x_FOR):0.1:max(x_FOR)])+par_DI(2,k),'-','LineWidth',2,'Color', map(4-k,:))
    %plot([min(x_FOR):0.1:max(x_FOR)],VPD_RESP(par_DI(:,k),[min(x_FOR):0.1:max(x_FOR)]),'-','LineWidth',2,'Color', map(4-k,:))
end
colormap(brewermap([10],'BrBg'))
ylim([0 150]);
xlim([0 2.3]);
cbr=colorbar;
set(cbr,'YTickLabel',(0.4:0.2:1).*100)
set(get(cbr,'ylabel'),'String', {'ET_m PET_m^{-1} [%]'},'FontSize',14);
axis square
ylabel('ET_m [mm]');
xlabel('VPD_m [kPa]');

subplot(1,2,2);
scatter(x_PTL,y_PTL,20,w_PTL,'filled')
caxis([0.4 1])
w_PTL(w_PTL>1)=1;
map=flipud(brewermap([4],'BrBg'));
for k=1:3;
    if k==1;
        ind=find(w_PTL<=0.5);
    else
        ind=find(w_PTL<=0.5+0.25*k-0.25 & w_PTL>0.5+0.25*k-0.5);
    end
    hold on
    n_PTL(k)=length(ind);
    for n = 1:1000;
        
        if n==1;
            %par_DI_B(:,n)= nlinfit(x_PTL(ind),y_PTL(ind),VPD_RESP,[0.7 130],opts);
            par_DI_B(:,n)= polyfit(x_PTL(ind),y_PTL(ind),1);
            unc_sl(:,n)=par_DI_B(1,n).*([min(x_PTL):0.1:max(x_PTL)])+par_DI_B(2,n);
            %unc_sl(:,n)=VPD_RESP(par_DI_B(:,n),[min(x_PTL):0.1:max(x_PTL)]);
        else
            x=x_PTL(ind);
            y=y_PTL(ind);
            index = randsample(1:1:length(x),length(x),'true');
            %par_DI_B(:,n)= nlinfit(x(index),y(index),VPD_RESP,[0.7 130],opts);
            par_DI_B(:,n)= polyfit(x(index),y(index),1);
            %unc_sl(:,n)=VPD_RESP(par_DI_B(:,n),[min(x_PTL):0.1:max(x_PTL)]);
            unc_sl(:,n)=par_DI_B(1,n).*([min(x_PTL):0.1:max(x_PTL)])+par_DI_B(2,n);
        end
    end
    par_DI_PTL(:,k)= par_DI_B(:,1);
    par_DI_CI_PTL{k}=prctile(par_DI_B',[2.5 97.5]);
    %par_DI_PTL(:,k)= nlinfit(x_PTL(ind),y_PTL(ind),VPD_RESP,[0.7 130],opts);
     %par_DI_PTL(:,k)= polyfit(x_PTL(ind),y_PTL(ind),1);
     plotshaded([min(x_PTL):0.1:max(x_PTL)],prctile(unc_sl',[2.5 97.5]),'r');
     plot([min(x_PTL):0.1:max(x_PTL)],par_DI_PTL(1,k).*([min(x_PTL):0.1:max(x_PTL)])+par_DI_PTL(2,k),'-','LineWidth',2,'Color', map(4-k,:))
    %plot([min(x_PTL):0.1:max(x_PTL)],VPD_RESP(par_DI(:,k),[min(x_PTL):0.1:max(x_PTL)]),'-','LineWidth',2,'Color', map(4-k,:))
end
colormap(brewermap([10],'BrBg'))
ylim([0 150]);
xlim([0 2.3]);
cbr=colorbar;
set(cbr,'YTickLabel',(0.4:0.2:1).*100)
set(get(cbr,'ylabel'),'String', {'ET_m PET_m^{-1} [%]'},'FontSize',14);
axis square
ylabel('ET_m [mm]');
xlabel('VPD_m [kPa]');
%%
figure,
for k=1:3;
 plot([min(x_PTL):0.1:max(x_PTL)],par_DI_PTL(1,k).*([min(x_PTL):0.1:max(x_PTL)])+par_DI_PTL(2,k),'-','LineWidth',2,'Color', map(4-k,:))
 hold on
 plot([min(x_FOR):0.1:max(x_FOR)],par_DI(1,k).*([min(x_FOR):0.1:max(x_FOR)])+par_DI(2,k),'--','LineWidth',2,'Color', map(4-k,:))
end
%% plot VPD-ET sensitivity for all sites with individual slopes
% derive slopes for annual dryness index
x_FOR=VPD_MNT_MID_FOR;
y_FOR=LE_MNT_FOR;
z_FOR=SITE_FOR;
w_FOR=RN_MNT_FOR;
t_FOR=MNT_FOR;
i=find(~isnan(x_FOR) & ~isnan(y_FOR) & (t_FOR>=5 & t_FOR<=9)); % & ~isnan(z_FOR));
%i=find(~isnan(x_FOR) & ~isnan(y_FOR) & ~isnan(z_FOR) & ~isnan(w_FOR));
x_FOR=x_FOR(i);
y_FOR=y_FOR(i);
z_FOR=z_FOR(i);
w_FOR=w_FOR(i);

x_PTL=VPD_MNT_MID_PTL;
y_PTL=LE_MNT_PTL;
z_PTL=SITE_PTL;
w_PTL=RN_MNT_PTL;
t_PTL=MNT_PTL;
i=find(~isnan(x_PTL) & ~isnan(y_PTL) & (t_PTL>=5 & t_PTL<=9));% & ~isnan(z_PTL));
%i=find(~isnan(x_PTL) & ~isnan(y_PTL) & ~isnan(z_PTL) & ~isnan(w_PTL));
x_PTL=x_PTL(i);
y_PTL=y_PTL(i);
z_PTL=z_PTL(i);
w_PTL=w_PTL(i);

figure,
subplot(1,2,1)
for k = 1:length(unique(z_PTL));
    %plot(x_PTL(z_PTL==k),y_PTL(z_PTL==k),'.','MarkerSize',20);
    hold on
    x=x_PTL(z_PTL==k);
    y=y_PTL(z_PTL==k);
    w=w_PTL(z_PTL==k);
    if length(y>3)
        
    %pr=polyfit(w,y,1);
    %y=y-(pr(1)*w+pr(2));
    if length(x)<=3;
        
        R2_PTL(k)=NaN;
        pval_PTL(k)=NaN;
        slope_PTL(k)=NaN;
    else
        pr_PTL(:,k) = nlinfit(x,y,VPD_RESP,[70 125],opts);
         [R p]=corrcoef(VPD_RESP(pr(:,k),x),y);
        %[R p]=corrcoef(x,y);
        R2_PTL(k)=R(1,2)^2;
        pval_PTL(k)=p(1,2);
        if p(1,2)<0.05;
            par=polyfit(x,y,1);
            %plot([min(x):0.1:max(x)],par(1).*[min(x):0.1:max(x)]+par(2),'-','LineWidth',2)
            hold on
            plot([min(x):0.1:max(x)],VPD_RESP(pr_PTL(:,k),[min(x):0.1:max(x)]),'-k','LineWidth',2)
            k
            %pause(3)
            %close all
            slope_PTL(k)=par(1);
        else
            slope_PTL(k)=NaN;
        end
    end
    end
end
axis square
%%
clear pr
subplot(1,2,1)
for k = 1:length(unique(z_FOR));
    %plot(x_FOR(z_FOR==k),y_FOR(z_FOR==k),'.','MarkerSize',20);
    hold on
    x=x_FOR(z_FOR==k);
    y=y_FOR(z_FOR==k);
    w=w_FOR(z_FOR==k);
    if length(y>3)
    %pr=polyfit(w,y,1);
    %y=y-(pr(1)*w+pr(2));
    if length(x)<=3;
        
        R2_FOR(k)=NaN;
        pval_FOR(k)=NaN;
        slope_FOR(k)=NaN;
    else
        pr(:,k) = nlinfit(x,y,VPD_RESP,[70 125],opts);
        [R p]=corrcoef(VPD_RESP(pr(:,k),x),y);
        %[R p]=corrcoef(x,y);
        R2_FOR(k)=R(1,2)^2;
        pval_FOR(k)=p(1,2);
        if p(1,2)<0.05;
            par=polyfit(x,y,1);
            %plot([min(x):0.1:max(x)],par(1).*[min(x):0.1:max(x)]+par(2),'-','LineWidth',2)
            hold on
            plot([min(x):0.1:max(x)],VPD_RESP(pr(:,k),[min(x):0.1:max(x)]),'-','LineWidth',2,'Color',[0.5 0.5 0.5])
            k
            axis square
           % pause(2)
            %close all
            slope_FOR(k)=par(1);
        else
            slope_FOR(k)=NaN;
        end
    end
    end
end


%%
x_PTL=x_PTL(z_PTL<=thr);
y_PTL=y_PTL(z_PTL<=thr);
x_FOR=x_FOR(z_FOR<=thr);
y_FOR=y_FOR(z_FOR<=thr);
R=corrcoef(x_FOR,y_FOR);
R(1,2)^2

R=corrcoef(x_PTL,y_PTL);
R(1,2)^2
for k=1:21;
    i=find(x_PTL>k./10 & x_PTL<=k./10+0.1);
    x_PTL_MN(k)=nanmean(x_PTL(i));
    y_PTL_MN(k)=nanmean(y_PTL(i));
    
    i=find(x_FOR>k./10 & x_FOR<=k./10+0.1);
    x_FOR_MN(k)=nanmean(x_FOR(i));
    y_FOR_MN(k)=nanmean(y_FOR(i));
end

 plot(x_FOR_MN,y_FOR_MN,'-');
 hold on

 plot(x_PTL_MN,y_PTL_MN,'-');
%%
x_FOR=VPD_MNT_MID_FOR;
y_FOR=LE_MNT_FOR;
z_FOR=P_PET_FOR';
i=find(~isnan(x_FOR) & ~isnan(y_FOR) & ~isnan(z_FOR));
x_FOR=x_FOR(i);
y_FOR=y_FOR(i);
z_FOR=z_FOR(i);

x_PTL=VPD_MNT_MID_PTL;
y_PTL=LE_MNT_PTL;
z_PTL=P_PET_PTL';
i=find(~isnan(x_PTL) & ~isnan(y_PTL) & ~isnan(z_PTL));
x_PTL=x_PTL(i);
y_PTL=y_PTL(i);
z_PTL=z_PTL(i);

x_PTL=x_PTL(z_PTL>thr);
y_PTL=y_PTL(z_PTL>thr);
x_FOR=x_FOR(z_FOR>thr);
y_FOR=y_FOR(z_FOR>thr);
R=corrcoef(x_FOR,y_FOR);
R(1,2)^2

R=corrcoef(x_PTL,y_PTL);
R(1,2)^2
for k=1:21;
    i=find(x_PTL>k./10 & x_PTL<=k./10+0.1);
    x_PTL_MN(k)=nanmean(x_PTL(i));
    y_PTL_MN(k)=nanmean(y_PTL(i));
    
    i=find(x_FOR>k./10 & x_FOR<=k./10+0.1);
    x_FOR_MN(k)=nanmean(x_FOR(i));
    y_FOR_MN(k)=nanmean(y_FOR(i));
end
subplot(1,2,2)
 plot(x_FOR_MN,y_FOR_MN,'-');
 hold on
 subplot(1,2,1)
 plot(x_PTL_MN,y_PTL_MN,'-');
%% read annual and growing season data
clearvars -except SITE_NAME
LE_ANN_PTL=[];
LE_ANN_FOR=[];

RN_GS_PTL=[];
RN_GS_FOR=[];

SWIN_GS_PTL=[];
SWIN_GS_FOR=[];

P_ANN_PTL=[];
P_ANN_FOR=[];

LE_GS_PTL=[];
LE_GS_FOR=[];
VPD_ANN_PTL=[];
VPD_ANN_FOR=[];

VPD_MID_PTL=[];
VPD_MID_FOR=[];

SITE_PTL=[];
SITE_FOR=[];
LAT_FOR=[];
ANN_FOR=[];
LON_FOR=[];

LAT_PTL=[];
ANN_PTL=[];
LON_PTL=[];
lng=length(SITE_NAME);


%%
%cd('/Users/manuelhelbig/Desktop/Mac_PDF_Boreal/ANN/');
cd('D:\Mac_PDF_Boreal\ANN\');
for k = 1:lng; %length(SITE_NAME);
    SITE_NAME{k}
    data=importdata(['ANN_' SITE_NAME{k} '.csv']);
    if nanmean(data(:,8))==1
%         subplot(1,2,1)
%         plot(data(:,4),data(:,2),'.','MarkerSize',24);
%         hold on
%         axis square
%         pause()
        LAT_PTL=vertcat(LAT_PTL,data(:,10));
        ANN_PTL=vertcat(ANN_PTL,data(:,1));
        LON_PTL=vertcat(LON_PTL,data(:,11));
        LE_ANN_PTL=vertcat(LE_ANN_PTL,data(:,2));
        SWIN_GS_PTL=vertcat(SWIN_GS_PTL,data(:,6));
        RN_GS_PTL=vertcat(RN_GS_PTL,data(:,5));
        P_ANN_PTL=vertcat(P_ANN_PTL,data(:,9));
        VPD_ANN_PTL=vertcat(VPD_ANN_PTL,data(:,4));
        if length(data(1,:))<12;
            VPD_MID_PTL=vertcat(VPD_MID_PTL,NaN(size(data(:,4))));
        else
            VPD_MID_PTL=vertcat(VPD_MID_PTL,data(:,12));
        end
        
        LE_GS_PTL=vertcat(LE_GS_PTL,data(:,3));
        SITE_PTL=vertcat(SITE_PTL,ones(size(data(:,3))).*k);
    else
%         subplot(1,2,2)
%         plot(data(:,4),data(:,2),'.','MarkerSize',24);
%         hold on
%         axis square
%         pause()
        LAT_FOR=vertcat(LAT_FOR,data(:,10));
        ANN_FOR=vertcat(ANN_FOR,data(:,1));
        LON_FOR=vertcat(LON_FOR,data(:,11));
        LE_ANN_FOR=vertcat(LE_ANN_FOR,data(:,2));
        P_ANN_FOR=vertcat(P_ANN_FOR,data(:,9));
        VPD_ANN_FOR=vertcat(VPD_ANN_FOR,data(:,4));
        if length(data(1,:))<12;
            VPD_MID_FOR=vertcat(VPD_MID_FOR,NaN(size(data(:,4))));
        else
            VPD_MID_FOR=vertcat(VPD_MID_FOR,data(:,12));
        end
        LE_GS_FOR=vertcat(LE_GS_FOR,data(:,3));
        SWIN_GS_FOR=vertcat(SWIN_GS_FOR,data(:,6));
        RN_GS_FOR=vertcat(RN_GS_FOR,data(:,5));
        SITE_FOR=vertcat(SITE_FOR,ones(size(data(:,3))).*k);
    end
end
%% get precip data for individual years
%cd('/Users/manuelhelbig/Desktop/CRU/');
cd('D:\CRU\');
path = '/Users/manuelhelbig/Desktop/CRU/';
path = 'D:\CRU\';
finfo = dir([path '*.nc']);
for k = 1:5;
    if k == 1;
    lat = ncread([path,finfo(k).name],'lat');
    lon = ncread([path,finfo(k).name],'lon');
    end
    if k <= 4;
        PRE(:,:,k*120-119:k*120) = ncread([path,finfo(k*4-2).name],'pre');
    else
        PRE(:,:,k*120-119:k*120-119+83) = ncread([path,finfo(k*4-2).name],'pre');
    end
    %{
    if k <= 4;
        TMP(:,:,k*120-119:k*120) = ncread([path,finfo(k*4-1).name],'tmp');
    else
        TMP(:,:,k*120-119:k*120-119+83) = ncread([path,finfo(k*4-1).name],'tmp');
    end
    %}
    if k <= 4;
        PET(:,:,k*120-119:k*120) = ncread([path,finfo(k*4-3).name],'pet');
    else
        PET(:,:,k*120-119:k*120-119+83) = ncread([path,finfo(k*4-3).name],'pet');
    end

end

MM=[];
YR=[];
for k=1:564/12;
    MM=horzcat(MM,1:12);
    YR=horzcat(YR,ones(1,12).*(1970+k));
end
MM_D=[31 28 31 30 31 30 31 31 30 31 30 31];
for k = 1:564;
    PET(:,:,k)=PET(:,:,k).*MM_D(MM(k));
end

for k = 1:564/12;
    PRE_ANN(:,:,k)=nansum(PRE(:,:,YR==1970+k),3);
    PRE_GS(:,:,k)=nansum(PRE(:,:,(YR==1970+k & (MM>=5 & MM<=9))),3);
    PET_ANN(:,:,k)=nansum(PET(:,:,YR==1970+k),3);
    PET_GS(:,:,k)=nansum(PET(:,:,(YR==1970+k & (MM>=5 & MM<=9))),3);
end
YY=1971:2017;
%%
for k=1:length(LON_FOR);
    lon_i=LON_FOR(k);
    lat_i=LAT_FOR(k);
    i1=find((lon-lon_i).^2==min((lon-lon_i).^2),1,'first');
    i2=find((lat-lat_i).^2==min((lat-lat_i).^2),1,'first');
    if isempty(find(YY==ANN_FOR(k)))
        PRE_FOR(k)=NaN;
        PRE_GS_FOR(k)=NaN;
        PET_FOR(k)=NaN;
        PET_GS_FOR(k)=NaN;
        
        PRE_FOR_CL(k)=NaN;
        PRE_GS_FOR_CL(k)=NaN;
        PET_FOR_CL(k)=NaN;
        PET_GS_FOR_CL(k)=NaN;
    else
        PRE_FOR(k)=PRE_ANN(i1,i2,find(YY==ANN_FOR(k)));
        PRE_GS_FOR(k)=PRE_GS(i1,i2,find(YY==ANN_FOR(k)));
        PET_FOR(k)=PET_ANN(i1,i2,find(YY==ANN_FOR(k)));
        PET_GS_FOR(k)=PET_GS(i1,i2,find(YY==ANN_FOR(k)));
        
        PRE_FOR_CL(k)=mean(PRE_ANN(i1,i2,find(YY>=2000 & YY<=2017)));
        PRE_GS_FOR_CL(k)=mean(PRE_GS(i1,i2,find(YY>=2000 & YY<=2017)));
        PET_FOR_CL(k)=mean(PET_ANN(i1,i2,find(YY>=2000 & YY<=2017)));
        PET_GS_FOR_CL(k)=mean(PET_GS(i1,i2,find(YY>=2000 & YY<=2017)));
    end
end

for k=1:length(LON_PTL);
    lon_i=LON_PTL(k);
    lat_i=LAT_PTL(k);
    i1=find((lon-lon_i).^2==min((lon-lon_i).^2),1,'first');
    i2=find((lat-lat_i).^2==min((lat-lat_i).^2),1,'first');
    if isempty(find(YY==ANN_PTL(k)))
        PRE_PTL(k)=NaN;
        PRE_GS_PTL(k)=NaN;
        PET_PTL(k)=NaN;
        PET_GS_PTL(k)=NaN;
        
        PRE_PTL_CL(k)=NaN;
        PRE_GS_PTL_CL(k)=NaN;
        PET_PTL_CL(k)=NaN;
        PET_GS_PTL_CL(k)=NaN;
    else
        PRE_PTL(k)=PRE_ANN(i1,i2,find(YY==ANN_PTL(k)));
        PRE_GS_PTL(k)=PRE_GS(i1,i2,find(YY==ANN_PTL(k)));
        PET_PTL(k)=PET_ANN(i1,i2,find(YY==ANN_PTL(k)));
        PET_GS_PTL(k)=PET_GS(i1,i2,find(YY==ANN_PTL(k)));
        
        PRE_PTL_CL(k)=mean(PRE_ANN(i1,i2,find(YY>=2000 & YY<=2017)));
        PRE_GS_PTL_CL(k)=mean(PRE_GS(i1,i2,find(YY>=2000 & YY<=2017)));
        PET_PTL_CL(k)=mean(PET_ANN(i1,i2,find(YY>=2000 & YY<=2017)));
        PET_GS_PTL_CL(k)=mean(PET_GS(i1,i2,find(YY>=2000 & YY<=2017)));
    end
end
%% get recent P-PET balance
PET_00_17=nanmean(PET_ANN(:,:,30:47),3);
P_00_17=nanmean(PRE_ANN(:,:,30:47),3);

%%
clear x_FOR y_FOR z_FOR x_PTL y_PTL z_PTL 
x_FOR=VPD_MID_FOR;
y_FOR=LE_GS_FOR;
z_FOR=PRE_GS_FOR;
l_FOR=RN_GS_FOR;
r_FOR=SWIN_GS_FOR;
%ind=find(~isnan(x_FOR)&~isnan(y_FOR)&~isnan(z_FOR)'&~isnan(l_FOR));
%ind=find(~isnan(x_FOR)&~isnan(y_FOR)&~isnan(l_FOR)&~isnan(z_FOR)');
ind=find(~isnan(y_FOR));
x_FOR=x_FOR(ind);
y_FOR=y_FOR(ind);
z_FOR=z_FOR(ind)';
l_FOR=l_FOR(ind);
r_FOR=r_FOR(ind);

x_PTL=VPD_MID_PTL;
y_PTL=LE_GS_PTL;
z_PTL=PRE_GS_PTL;
l_PTL=RN_GS_PTL;
r_PTL=SWIN_GS_PTL;
%ind=find(~isnan(x_PTL)&~isnan(y_PTL)&~isnan(z_PTL)'&~isnan(l_PTL));
ind=find(~isnan(y_PTL));
%ind=find(~isnan(x_PTL)&~isnan(y_PTL)&~isnan(l_PTL)&~isnan(z_PTL)');
x_PTL=x_PTL(ind);
y_PTL=y_PTL(ind);
z_PTL=z_PTL(ind)';
l_PTL=l_PTL(ind);
r_PTL=r_PTL(ind);

[pval, observeddifference, effectsize] = permutationTest(y_FOR,y_PTL,999);

%%

hist(x_PTL);
hold on
hist(x_FOR);
%%
thr=1;
for k = 1:2;
    if k==1;
        ind=find(z_PTL>thr);
    else
        ind=find(z_PTL<=thr);
    end
    
x=x_PTL(ind);
y=y_PTL(ind);
z=l_PTL(ind);

par_Rn = polyfit(z,y,1);
y=y-(par_Rn(1).*z+par_Rn(2));

ET_P_PTL(k)=nanmean(y);
ET_P_STD_PTL(k)=nanstd(y);
clear par unc_sl
    for n = 1:1000;
        
        if n==1;
            par(:,n) = polyfit(x,y,1);
            unc_sl(:,n)=par(1).*(0.2:0.1:1)+par(2);
        else
            index = randsample(1:1:length(x),length(x),'true');
            par(:,n) = polyfit(x(index),y(index),1);
            unc_sl(:,n)=par(1,n).*(0.2:0.1:1)+par(2,n);
        end
    end
    sl_PTL(k)=par(1,1);
    low_CI_PTL(k)=prctile(par(1,:),2.5);
    up_CI_PTL(k)=prctile(par(1,:),97.5);
    std_PTL(k)=std(par(1,:));
    %plot(par(1,1).*x+par(2,1),y-(par(1,1).*x+par(2,1)),'o')
    %pause()
    if k==1;
        par_WET_PTL=par;
    else
        par_DRY_PTL=par;
    end
end
%%
for k = 1:2;
    if k==1;
        ind=find(z_FOR>thr);
    else
        ind=find(z_FOR<=thr);
    end
x=x_FOR(ind);
y=y_FOR(ind);
z=l_FOR(ind);

par_Rn = polyfit(z,y,1);
%par_Rn
y=y-(par_Rn(1).*z+par_Rn(2));

ET_P_FOR(k)=nanmean(y);
ET_P_STD_FOR(k)=nanstd(y);
clear par unc_sl
    for n = 1:1000;
        
        if n==1;
            par(:,n) = polyfit(x,y,1);
            unc_sl(:,n)=par(1).*(0.2:0.1:1)+par(2);
        else
            index = randsample(1:1:length(x),length(x),'true');
            par(:,n) = polyfit(x(index),y(index),1);
            unc_sl(:,n)=par(1,n).*(0.2:0.1:1)+par(2,n);
        end
    end
    sl_FOR(k)=par(1,1);
    low_CI_FOR(k)=prctile(par(1,:),2.5);
    up_CI_FOR(k)=prctile(par(1,:),97.5);
    std_FOR(k)=std(par(1,:));
    %plot(par(1,1).*x+par(2,1),y-(par(1,1).*x+par(2,1)),'o')
    %pause()
     if k==1;
        par_WET_FOR=par;
    else
        par_DRY_FOR=par;
    end

end
%% create alternative figure with climatological P-PET
% temporal vs spatial relationship
figure,
subplot(1,2,1);
sl_PTL=par_PTL{1}(:,1);
l1=bar(0.25,sl_PTL(1),0.3,'r')      
hold on

er = errorbar(0.25,sl_PTL(1),sl_PTL(1)-prctile(par_PTL{1}(1,:),2.5),prctile(par_PTL{1}(1,:),97.5)-sl_PTL(1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

hold on

sl_FOR=par_FOR{1}(:,1);
l2=bar(0.5,sl_FOR(1),0.3,'r')                
hold on

er = errorbar(0.5,sl_FOR(1),sl_FOR(1)-prctile(par_FOR{1}(1,:),2.5),prctile(par_FOR{1}(1,:),97.5)-sl_FOR(1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
legend([l1 l2],{'PTL','FOR'});
ylabel('b [??]');
ax = gca;
set(gca,'FontSize',16);

subplot(1,2,2);
hold on
bar(0.85,sl_PTL(2),0.3,'b')                
hold on

er = errorbar(0.85,sl_PTL(2),sl_PTL(2)-prctile(par_PTL{1}(2,:),2.5),prctile(par_PTL{1}(2,:),97.5)-sl_PTL(2));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off



hold off
hold on
bar(1.1,sl_FOR(2),0.3,'b')                
hold on

er = errorbar(1.1,sl_FOR(2),sl_FOR(2)-prctile(par_FOR{1}(2,:),2.5),prctile(par_FOR{1}(2,:),97.5)-sl_FOR(2));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  


hold on
ylabel('ET_{ref} [mm]');
%plot([1.3 1.3],[-200 600],'--k','LineWidth',2);
ax = gca;
set(gca,'FontSize',16);
%%
%{
load('/Volumes/Seagate Backup Plus Drive/CMIP5/ET-VPD/STATS_ESM.mat');
l1=bar(1.45:0.75:6.7,par_WET(1,:),0.3,'r');
hold on
l1=bar(1.7:0.75:6.95,par_DRY(1,:),0.3,'b');
xlim([0 7.5]);
ax.XTick = [0.3:0.725:7.05];
%xlim([0 1.3]);
set(gca,'XTickLabel',{'PTL','FOR',ESM_NAME{1},ESM_NAME{2},ESM_NAME{3}...
    ,ESM_NAME{4},ESM_NAME{5},ESM_NAME{6},ESM_NAME{7},ESM_NAME{8}});
%}
%%
%{
low_edg=250:25:800;
up_edg=550:25:1100;
clear sl_PTL low_CI_PTL up_CI_PTL sl_FOR low_CI_FOR up_CI_FOR
for k = 1:length(low_edg);
    if k<length(low_edg) & k~=1;
        ind=find(z_PTL>low_edg(k) & z_PTL<=up_edg(k));
        ind_FOR=find(z_FOR>low_edg(k) & z_FOR<=up_edg(k));
    elseif k==1;
        ind=find(z_PTL<=up_edg(k));
        ind_FOR=find(z_FOR<=up_edg(k));
    else
        ind=find(z_PTL>low_edg(k));
        ind_FOR=find(z_FOR>low_edg(k));
    end
    lPTL(k)=length(y_PTL(ind));
    lFOR(k)=length(y_FOR(ind_FOR));
    x=x_PTL(ind);
    y=y_PTL(ind);
    ET_P_PTL(k)=nanmean(y);
    ET_P_STD_PTL(k)=nanstd(y);
    clear par unc_sl
    for n = 1:1000;
        
        if n==1;
            par(:,n) = polyfit(x,y,1);
            unc_sl(:,n)=par(1).*(0.2:0.1:1)+par(2);
        else
            index = randsample(1:1:length(x),length(x),'true');
            par(:,n) = polyfit(x(index),y(index),1);
            unc_sl(:,n)=par(1,n).*(0.2:0.1:1)+par(2,n);
        end
    end
    sl_PTL(k)=par(1,1);
    low_CI_PTL(k)=prctile(par(1,:),2.5);
    up_CI_PTL(k)=prctile(par(1,:),97.5);
    std_PTL(k)=std(par(1,:));
    
    x=x_FOR(ind_FOR);
    y=y_FOR(ind_FOR);
    
    ET_P_FOR(k)=nanmean(y);
    ET_P_STD_FOR(k)=nanstd(y);
    clear par unc_sl
    for n = 1:1000;
        
        if n==1;
            par(:,n) = polyfit(x,y,1);
            unc_sl(:,n)=par(1).*(0.2:0.1:1)+par(2);
        else
            index = randsample(1:1:length(x),length(x),'true');
            par(:,n) = polyfit(x(index),y(index),1);
            unc_sl(:,n)=par(1,n).*(0.2:0.1:1)+par(2,n);
        end
    end
    sl_FOR(k)=par(1,1);
    low_CI_FOR(k)=prctile(par(1,:),2.5);
    up_CI_FOR(k)=prctile(par(1,:),97.5);
    std_FOR(k)=std(par(1,:));
%     plot(x,y,'o')
%     hold on
%     plotshaded((0.2:0.1:1),prctile(unc_sl',[2.5 97.5]),'r');
%     plot((0.2:0.1:1),unc_sl(:,1),'-');
%     pause()
%     close all
end
%%
figure,
subplot(1,2,1)
plotshaded(375:25:925,[(sl_PTL-std_PTL)' (sl_PTL+std_PTL)'],'r');
%plotshaded(375:25:925,[low_CI_PTL' up_CI_PTL'],'r');
hold on
plot(375:25:925,sl_PTL,'-');

plotshaded(375:25:925,[(sl_FOR-std_FOR)' (sl_FOR+std_FOR)'],'r');
%plotshaded(375:25:925,[low_CI_FOR' up_CI_FOR'],'r');
hold on
plot(375:25:925,sl_FOR,'-');


xlim([350 950]);
axis square
refline(0,0);
ylim([-75 625]);

subplot(1,2,2)
plotshaded(375:25:925,[(ET_P_PTL-ET_P_STD_PTL)' (ET_P_PTL+ET_P_STD_PTL)'],'r');
hold on
plot(375:25:925,ET_P_PTL,'-');

plotshaded(375:25:925,[(ET_P_FOR-ET_P_STD_FOR)' (ET_P_FOR+ET_P_STD_FOR)'],'r');
hold on
plot(375:25:925,ET_P_FOR,'-');
axis square
xlim([350 950]);
refline(0,0);
ylim([200 525]);
%}
%%
clear x_FOR y_FOR z_FOR x_PTL y_PTL z_PTL 
x_FOR=VPD_MID_FOR;
y_FOR=LE_GS_FOR;
z_FOR=PRE_FOR-PET_FOR;
z_FOR1=PRE_FOR./PET_FOR;
ind=find(~isnan(x_FOR)&~isnan(y_FOR)&~isnan(z_FOR'));
x_FOR=x_FOR(ind);
y_FOR=y_FOR(ind);
z_FOR=z_FOR(ind)';
z_FOR1=z_FOR1(ind)';

x_PTL=VPD_MID_PTL;
y_PTL=LE_GS_PTL;
z_PTL=PRE_PTL-PET_PTL;
z_PTL1=PRE_PTL./PET_PTL;
ind=find(~isnan(x_PTL)&~isnan(y_PTL)&~isnan(z_PTL'));
x_PTL=x_PTL(ind);
y_PTL=y_PTL(ind);
z_PTL=z_PTL(ind)';
z_PTL1=z_PTL1(ind)';

[pval, observeddifference, effectsize] = permutationTest(y_FOR,y_PTL,999);
%
figure,
subplot(1,2,1);

[R pval] = corrcoef(x_PTL(z_PTL>0),y_PTL(z_PTL>0));
R(1,2)^2
pval(1,2)

clear par unc_sl
x=x_PTL(z_PTL1>thr);
y=y_PTL(z_PTL1>thr);
for k = 1:1000;
    if k==1;
        par(:,k) = polyfit(x,y,1);
        unc_sl(:,k)=par(1).*(0.2:0.1:1.85)+par(2);
    else
        index = randsample(1:1:length(x),length(x),'true');
        par(:,k) = polyfit(x(index),y(index),1);
        unc_sl(:,k)=par(1,k).*(0.2:0.1:1.85)+par(2,k);
    end
end
par = polyfit(x,y,1);

%scatter(x_PTL,y_PTL,z_PTL./8);
z_PTL(z_PTL<=-150)=-225;
z_PTL(z_PTL>-150 & z_PTL<=0)=-75;
z_PTL(z_PTL>0 & z_PTL<=150)=75;
z_PTL(z_PTL>150 & z_PTL<=300)=225;
z_PTL(z_PTL>300)=375;
z_PTL=z_PTL+525;
sze=[300 450 600 750 900];
bubsizes = unique(round((z_PTL./5)./10).*10)';
legentry=cell(size(bubsizes));
for ind = 1:numel(bubsizes)
   bubleg(ind) = plot(0,0,'bo','markersize',sqrt(bubsizes(ind)),'MarkerFaceColor','blue');
   hold on
   set(bubleg(ind),'visible','off')
   legentry{ind} = num2str(sze(ind));
end

h = scatter(x_PTL(z_PTL<75+525),y_PTL(z_PTL<75+525),z_PTL(z_PTL<75+525)./5,'r','MarkerFaceColor','red');
hold on

scatter(x_PTL(z_PTL>=75+525),y_PTL(z_PTL>=75+525),z_PTL(z_PTL>=75+525)./5,'b','MarkerFaceColor','blue');
%legend(legentry)

plotshaded(0.2:0.1:1.85,prctile(unc_sl',[2.5 97.5]),'r');
hold on
plot(0.2:0.1:1.85,par(1).*(0.2:0.1:1.85)+par(2),'-','LineWidth',2);
hold on
%xlabel('Mean VPD_{gs} [kPa]');
ylabel('Annual ET [mm]');
text(0.4,120,['P-PET > 0: R^2 = ' num2str(R(1,2)^2) ', p < 0.001, n = ' num2str(length(x))]);

text(0.4,100,['P-PET < 0: R^2 = n.s., n = ' num2str(length(x_PTL(z_PTL<75+525)))]);
axis square
ylim([100 650]);
xlim([0.2 1.75]);
title('PTL');
set(gca,'FontSize',16);
clear x y
clear pval R
subplot(1,2,2);
%plot(x_FOR,y_FOR,'o');
%scatter(x_FOR,y_FOR,z_FOR./8);
z_FOR(z_FOR<=-150)=-225;
z_FOR(z_FOR>-150 & z_FOR<=0)=-75;
z_FOR(z_FOR>0 & z_FOR<=150)=75;
z_FOR(z_FOR>150 & z_FOR<=300)=225;
z_FOR(z_FOR>300)=375;
z_FOR=z_FOR+525;
bubsizes = unique(round((z_FOR./5)./10).*10)';
legentry=cell(size(bubsizes));
for ind = 1:numel(bubsizes)
   bubleg(ind) = plot(0,0,'bo','markersize',sqrt(bubsizes(ind)),'MarkerFaceColor','blue');
   hold on
   set(bubleg(ind),'visible','off')
   legentry{ind} = num2str(sze(ind)-525);
end
h = scatter(x_FOR(z_FOR<75+525),y_FOR(z_FOR<75+525),z_FOR(z_FOR<75+525)./5,'r','MarkerFaceColor','red');
hold on
scatter(x_FOR(z_FOR>=75+525),y_FOR(z_FOR>=75+525),z_FOR(z_FOR>=75+525)./5,'b','MarkerFaceColor','blue');
lgd=legend(legentry)
title(lgd,'P_{ann}-PET_{ann} [mm]');
hold on

xlabel('Mean VPD_{max, gs} [kPa]');

title('FOR');
ylim([100 650]);
xlim([0.2 1.75]);
axis square
set(gca,'FontSize',16);

clear x_FOR y_FOR z_FOR
x_FOR=VPD_MID_FOR;
y_FOR=LE_GS_FOR;
z_FOR=PRE_FOR-PET_FOR;
z_FOR1=PRE_FOR./PET_FOR;
ind=find(~isnan(x_FOR)&~isnan(y_FOR)&~isnan(z_FOR'));
x_FOR=x_FOR(ind);
y_FOR=y_FOR(ind);
z_FOR=z_FOR(ind)';
z_FOR1=z_FOR1(ind)';

clear par unc_sl
x=x_FOR(z_FOR1<=thr);
y=y_FOR(z_FOR1<=thr);
for k = 1:1000;
    if k==1;
        par(:,k) = polyfit(x,y,1);
        unc_sl(:,k)=par(1).*(0.2:0.1:1.85)+par(2);
    else
        index = randsample(1:1:length(x),length(x),'true');
        par(:,k) = polyfit(x(index),y(index),1);
        unc_sl(:,k)=par(1,k).*(0.2:0.1:1.85)+par(2,k);
    end
end
[R pval] = corrcoef(x,y);
R(1,2)^2
pval(1,2)
text(0.4,100,['P-PET < 0: R^2 = ' num2str(R(1,2)^2) ', p = ' num2str(pval(1,2)) ', n = ' num2str(length(x))]);
hold on
plotshaded(0.2:0.1:1.85,prctile(unc_sl',[2.5 97.5]),'r');
hold on
plot(0.2:0.1:1.85,par(1,1).*(0.2:0.1:1.85)+par(2,1),'-','LineWidth',2);
hold on
clear par unc_sl

x=x_FOR(z_FOR1>thr);
y=y_FOR(z_FOR1>thr);
for k = 1:1000;
    if k==1;
        par(:,k) = polyfit(x,y,1);
        unc_sl(:,k)=par(1).*(0.2:0.1:1.85)+par(2);
    else
        index = randsample(1:1:length(x),length(x),'true');
        par(:,k) = polyfit(x(index),y(index),1);
        unc_sl(:,k)=par(1,k).*(0.2:0.1:1.85)+par(2,k);
    end
end
[R pval] = corrcoef(x,y);
R(1,2)^2
pval(1,2)
text(0.4,120,['P-PET > 0: R^2 = ' num2str(R(1,2)^2) ', p = ' num2str(pval(1,2)) ', n = ' num2str(length(x))]);
hold on
plotshaded(0.2:0.1:1.85,prctile(unc_sl',[2.5 97.5]),'r');
hold on
plot(0.2:0.1:1.85,par(1,1).*(0.2:0.1:1.85)+par(2,1),'-','LineWidth',2);
hold on
clear par unc_sl
%%
x_FOR=VPD_ANN_FOR;
y_FOR=LE_ANN_FOR;
z_FOR=PRE_FOR;
ind=find(~isnan(z_FOR')&~isnan(y_FOR));
x_FOR=x_FOR(ind);
y_FOR=y_FOR(ind);
z_FOR=z_FOR(ind)';

x_PTL=VPD_ANN_PTL;
y_PTL=LE_ANN_PTL;
z_PTL=PRE_PTL;
ind=find(~isnan(z_PTL')&~isnan(y_PTL));
x_PTL=x_PTL(ind);
y_PTL=y_PTL(ind);
z_PTL=z_PTL(ind)';

%%
figure,
subplot(1,2,1);

[R pval] = corrcoef(z_PTL,y_PTL);
R(1,2)^2
pval(1,2)
clear par unc_sl
for k = 1:1000;
    if k==1;
        par(:,k) = polyfit(z_PTL,y_PTL,1);
        unc_sl(:,k)=par(1).*(100:1000)+par(2);
    else
        index = randsample(1:1:length(z_PTL),length(z_PTL),'true');
        par(:,k) = polyfit(z_PTL(index),y_PTL(index),1);
        unc_sl(:,k)=par(1,k).*(100:1000)+par(2,k);
    end
end
plotshaded(100:1000,prctile(unc_sl',[2.5 97.5]),'r');
hold on
clear par
par = polyfit(z_PTL,y_PTL,1);
plotshaded(100:1000,prctile(unc_sl',[2.5 97.5]),'r');
hold on
plot(100:1000,par(1).*(100:1000)+par(2),'-','LineWidth',2);
%scatter(z_PTL,y_PTL,18,x_PTL,'filled');
plot(z_PTL,y_PTL,'o');
%plot(100:1000,100:1000,'-');
hold on
xlabel('Annual P [m]');
ylabel('Annual ET [mm]');
text(150,120,['R^2 = ' num2str(R(1,2)^2) ', p < 0.001, n = ' num2str(length(z_PTL))]);
axis square
ylim([100 650]);
xlim([100 1000]);
title('PTL');
set(gca,'FontSize',16);

subplot(1,2,2);
[R pval] = corrcoef(z_FOR,y_FOR);
R(1,2)^2
pval(1,2)
plot(z_FOR,y_FOR,'o');
%scatter(x_FOR,y_FOR,18,z_FOR,'filled');
hold on
clear par unc_sl
for k = 1:1000;
    if k==1;
        par(:,k) = polyfit(z_FOR,y_FOR,1);
        unc_sl(:,k)=par(1).*(100:1000)+par(2);
    else
        index = randsample(1:1:length(z_FOR),length(z_FOR),'true');
        par(:,k) = polyfit(z_FOR(index),y_FOR(index),1);
        unc_sl(:,k)=par(1,k).*(100:1000)+par(2,k);
    end
end
plotshaded(100:1000,prctile(unc_sl',[2.5 97.5]),'r');
hold on
clear par
par = polyfit(z_FOR,y_FOR,1);
plot(100:1000,par(1).*(100:1000)+par(2),'-','LineWidth',2);

text(150,120,['R^2 = ' num2str(R(1,2)^2) ', p < 0.001, n = ' num2str(length(z_FOR))]);
title('FOR');
ylim([100 650]);
xlim([100 1000]);
axis square
set(gca,'FontSize',16);