%% Preprocessing: creazione matrici "distanze", "dimensioni", "attivazioni"

load('WS_1_001.mat')
dimensioni=ssirl(1:68,:);
attivazioni=Mfret(1:68,:);
distanze=distance(1:68,:);
load('WS_1_002.mat')
dimensioni=cat(2,dimensioni,ssirl(1:68,:));
attivazioni=cat(2,attivazioni,Mfret(1:68,:));
distanze=cat(2,distanze,distance(1:68,:));
load('WS_1_003.mat')
dimensioni=cat(2,dimensioni,ssirl(1:68,:));
attivazioni=cat(2,attivazioni,Mfret(1:68,:));
distanze=cat(2,distanze,distance(1:68,:));
load('WS_2_001.mat')
dimensioni=cat(2,dimensioni,ssirl(1:68,:));
attivazioni=cat(2,attivazioni,Mfret(1:68,:)); % decidere se metterlo o no
distanze=cat(2,distanze,distance(1:68,:));
load('WS_2_002.mat')
dimensioni=cat(2,dimensioni,ssirl(1:68,:));
attivazioni=cat(2,attivazioni,Mfret(1:68,:));
distanze=cat(2,distanze,distance(1:68,:));
load('WS_2_003.mat')
dimensioni=cat(2,dimensioni,ssirl(1:68,:));
attivazioni=cat(2,attivazioni,Mfret(1:68,:));
distanze=cat(2,distanze,distance(1:68,:));
load('WS_2_004.mat')
dimensioni=cat(2,dimensioni,ssirl(1:68,:));
attivazioni=cat(2,attivazioni,Mfret(1:68,:));
distanze=cat(2,distanze,distance(1:68,:));
load('WS_3_002.mat')
dimensioni=cat(2,dimensioni,ssirl(1:68,:));
attivazioni=cat(2,attivazioni,Mfret(1:68,:));
distanze=cat(2,distanze,distance(1:68,:));
% load('WS_3_003.mat')
% dimensioni=cat(2,dimensioni,ssirl(1:68,:));
% attivazioni=cat(2,attivazioni,Mfret(1:68,:));
% distanze=cat(2,distanze,distance(1:68,:));
load('WS_3_004.mat')
dimensioni=cat(2,dimensioni,ssirl(1:68,:));
attivazioni=cat(2,attivazioni,Mfret(1:68,:));
distanze=cat(2,distanze,distance(1:68,:));
% load('WS_2_sonda_dmso+trf_Series_4.mat')
% dimensioni=cat(2,dimensioni,ssirl(1:68,:));
% attivazioni=cat(2,attivazioni,Mfret(1:68,:));
% distanze=cat(2,distanze,distance(1:68,:));
% load('WS_sonda_dmso+trf da sottrarre il pallozzo_Series_1.mat')
% dimensioni=cat(2,dimensioni,ssirl(1:68,:));
% attivazioni=cat(2,attivazioni,Mfret(1:68,:));
% distanze=cat(2,distanze,distance(1:68,:));
% load('WS_sonda_dmso+trf da sottrarre il pallozzo_Series_2.mat')
% dimensioni=cat(2,dimensioni,ssirl(1:68,:));
% attivazioni=cat(2,attivazioni,Mfret(1:68,:));
% distanze=cat(2,distanze,distance(1:68,:));

%% Binning in funzione della distanza dal nucleo

% figure,
% 
% for c=0:7
%     p1=c*50;
%     p2=c*50+50;
%     
%     clear dimensioneP dimensione dimensioneG attivazione attivazioneP attivazioneG attivazionePnan attivazioneGnan dimensionePnan dimensioneGnan
%     
%     piccole=0; % il nome piccole non ha significato, sono solo pigro
%     for i=1:5
%         p=find(distanze(i,:)<p2 & distanze(i,:)>p1);
%         if isempty(p)==0
%             for j=1:size(p,2)
%                 piccole(j,i)=p(1,j);
%             end
%         end
%     end
%     distanzeIni=0;
%     attivazioneD=0;
%     for i=1:size(piccole,2)
%         for k=1:nnz(piccole(:,i))
%             distanzeIni(k,i)=distanze(i,piccole(k,i));
%             attivazioneD(k,i)=attivazioni(i,piccole(k,i));
%         end
%     end
%     
%     distanzeIninan=distanzeIni;
%     attivazioneDnan=attivazioneD;
%     distanzeIninan(attivazioneDnan==0)=NaN;
%     attivazioneDnan(attivazioneDnan==0)=NaN;
%     
%     
%     Mdis=nanmean(distanzeIninan(:));
%     MattD=nanmean(attivazioneDnan(:));
%     
%     SDdis=nanstd(distanzeIninan(:));
%     SDattD=nanstd(attivazioneDnan(:));
%     
%     for i=1:size(distanzeIninan,2)
%         scatter(distanzeIninan(:,i)*62.5,attivazioneDnan(:,i),'b.'), hold on,
%     end
%     scatter(Mdis*62.5,MattD,'go'), hold on,
%     errorbar(Mdis*62.5,MattD,SDattD),
%     herrorbar(Mdis*62.5,MattD,SDdis*62.5)
%     
%     mediaDist(c+1)=Mdis;
%     mediaAtt(c+1)=MattD;
%     stdDist(c+1)=SDdis;
%     stdAtt(c+1)=SDattD;
%     numero(c+1)=nnz(distanzeIni);
% end
% 
% for i=1:8
% ser(i)=stdAtt(i)/sqrt(numero(i));
% end
% 
% binCent=[25,75,125,175,225,275,325,375];
% figure, bar(binCent*62.5,mediaAtt), hold on,
% errorbar(binCent*62.5,mediaAtt,ser)


for i=1:68
    for c=1:2
%         p1=c*50;
%         p2=c*50+50;
        if c==1
            p1=30;
            p2=100;
        end
        if c==2
            p1=100;
            p2=400;
        end
        p=0;
        p=find(distanze(i,:)<p2 & distanze(i,:)>p1);
        binned=0;
        if isempty(p)==0
            for j=1:size(p,2)
                binned(j,i)=p(1,j);  % non serve scrivere binned come matrice, tanto ad ogni bin viene sovrascritto (!?! e quelli che non vengono sovrascritti?)
            end
        end
        distBin=0;
        dimBinPix=0;
        attBin=0;
        dimBin=0;
        for ii=1:size(binned,2)            % funziona ma fa schifo, le "i" non servono
            for k=1:nnz(binned(:,ii))
                distBin(k,i)=distanze(i,binned(k,ii))*62.5/1000;
                attBin(k,i)=attivazioni(i,binned(k,ii));
                dimBin(k,i)=sqrt(dimensioni(i,binned(k,ii))/pi)*62.5;
                dimBinPix(k,i)=dimensioni(i,binned(k,ii));
            end
        end
        
        distBinNan=distBin;
        attBinNan=attBin;
        dimBinNan=dimBin;
%         dimBinPixNan=dimBinPix;
        distBinNan(attBinNan==0)=NaN;
        attBinNan(attBinNan==0)=NaN;
        dimBinNan(dimBinNan==0)=NaN;
%         dimBinPixNan(dimBinPixNan==0)=NaN;
        
        Mdis(i,c+1)=nanmean(distBinNan(:));
        Matt(i,c+1)=nanmean(attBinNan(:));
        Mdim(i,c+1)=nanmean(dimBinNan(:));
        MdimN(i,c+1)=nanmean(dimBinNan(:))/(sum(dimensioni(i,:))/nnz(dimensioni(i,:)));
        MdimNN(i,c+1)=nanmean(dimBinNan(:))/sum(dimensioni(i,:));
        
        SumDimPix(i,c+1)=sum(dimBinPix(:,i));
        
        tgMediaGlob(i)=sum(dimensioni(i,:))/nnz(dimensioni(i,:));
        tgGlob(i)=sum(dimensioni(i,:));
        
        SDdis(i,c+1)=nanstd(distBinNan(:));
        SDatt(i,c+1)=nanstd(attBinNan(:));
        SDdim(i,c+1)=nanstd(dimBinNan(:));
        
        numero(i,c+1)=nnz(distBin);
    end
end

for j=1:68
    for i=2:3
        serDim(j,i)=SDdim(j,i)/sqrt(numero(j,i));
        serAtt(j,i)=SDatt(j,i)/sqrt(numero(j,i));
    end
end

binCent=[25,75,125,175,225,275,325,375];
figure, bar(binCent*62.5,Matt(68,:)), hold on,
errorbar(binCent*62.5,Matt(68,:),serAtt(68,:))
figure, bar(binCent*62.5,Mdim(68,:)), hold on,
errorbar(binCent*62.5,Mdim(68,:),serDim(68,:))
hold on,
bar(binCent*62.5,Mdim(3,:)), hold on,
errorbar(binCent*62.5,Mdim(3,:),serDim(3,:))
mean(Mdim(2:6,:),1)
mean(Mdim(64:68,:),1)

figure, 
for a=1:8
    plot(Mdim(:,a)), hold on, errorbar(Mdim(:,a),serDim(:,a)), hold on,
end
figure, 
for a=1:8
    plot(Matt(:,a)), hold on,
end

figure,
plot(Mdim(2:68,2)/mean(Mdim(2:4,2)),'r', 'LineWidth', 2), hold on, errorbar(Mdim(2:68,2)/mean(Mdim(2:4,2)),serDim(2:68,2)/mean(Mdim(2:4,2))), hold on,
plot(Mdim(2:68,3)/mean(Mdim(2:4,3)),'b', 'LineWidth', 2), hold on, errorbar(Mdim(2:68,3)/mean(Mdim(2:4,3)),serDim(2:68,3)/mean(Mdim(2:4,3))),

figure,
plot(Matt(2:68,2),'r', 'LineWidth', 2), hold on, errorbar(Matt(2:68,2),serAtt(2:68,2)), hold on,
plot(Matt(2:68,3)+0.05,'b', 'LineWidth', 2), hold on, errorbar(Matt(2:68,3),serAtt(2:68,3)),


figure,
plot(Matt(2:68,2),'r', 'LineWidth', 2), hold on,
plot(Matt(2:68,3),'b', 'LineWidth', 2),

figure,
plot(Mdim(2:68,2),'r', 'LineWidth', 2), hold on,
plot(Mdim(2:68,3),'b', 'LineWidth', 2),

figure,
plot(MdimN(2:68,2),'r', 'LineWidth', 2), hold on,
plot(MdimN(2:68,3),'b', 'LineWidth', 2),

figure,
plot(MdimNN(2:68,2),'r', 'LineWidth', 2), hold on,
plot(MdimNN(2:68,3),'b', 'LineWidth', 2),

figure, plot(tgMediaGlob)
figure, plot(tgGlob)

% dimensione nei due bin normalizzata sul numero di pixel totali delle
% vescicole
dimTot=sum(SumDimPix,2); % pixel totali occupati da vescicole nei vari frame
for i=1:68
    SumDimPixN1(i)=SumDimPix(i,2)/dimTot(i);
    SumDimPixN2(i)=SumDimPix(i,3)/dimTot(i);
end
figure, plot(SumDimPixN1)
figure, plot(SumDimPixN2)
figure, plot(dimTot)

% dimensione nei due bin normalizzata sul numero totale di vescicole
numTot=sum(numero,2); % pixel totali occupati da vescicole nei vari frame
for i=1:68
    MdimN1(i)=Mdim(i,2)/numTot(i);
    MdimN2(i)=Mdim(i,3)/numTot(i);
end
figure, plot(MdimN1)
figure, plot(MdimN2)
for i=1:68
    NumN1(i)=numero(i,2)/numTot(i);
    NumN2(i)=numero(i,3)/numTot(i);
end
figure, plot(NumN1)
figure, plot(NumN2)
figure, plot(numTot)