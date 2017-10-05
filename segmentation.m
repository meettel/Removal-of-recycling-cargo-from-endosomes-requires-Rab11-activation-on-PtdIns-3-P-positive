clearvars -except RvettoreFrame Rprova
ftoty=0;
ftotc=0;
ftotr=0;
foy=0;
foc=0;
fore=0;
foyM=0;
focM=0;
foreM=0;

% PARAMETRI DA REGOLARE MANUALMENTE:

% NOTA: devi valutare manualmente anche la necessit� e il numero di filtri
% mediani e dilatazioni (righe 148 149 153 154 166 470 471 475 476 488) e settare la normalizzazione
% della fluorescenza (righe 287 288 289 687 688 689)

nXpos=260; % coordinate del punto dal quale si calcola la distanza delle vescicole
nYpos=80;

% r1x1=374; % punti che definiscono i 4 rettangoli scelti per calcolare le intensita` nel citoplasma
% r1x2=458;
% r1y1=42;
% r1y2=113;
% 
% r2x1=443;
% r2x2=510;
% r2y1=176;
% r2y2=265;
% 
% r3x1=342;
% r3x2=426;
% r3y1=419;
% r3y2=504;
% 
% r4x1=2;
% r4x2=73;
% r4y1=337;
% r4y2=405;

st=5; % numero stacks (4 o 5 nel caso siano di pi� o di meno ci sono da fare piccole modifiche)
minC=35; 
minR=35; % sono i valori minimi di intensit� del pixel mecessari per considerare segnale reale il dato (vedi punto di utilizzo per maggiori info)
gaussD=[3 3]; % dimensione del filtro gaussiano
gaussS=0.95; % deviazione standard del filtro gaussiano
thre1=0.52; % soglia del thresholding che trasforma l'immagine da segmentare in binaria
thre2=0.78;
contr=0.85; % upper limit contrasto
minDim=60; % dimensione minima in pixel che devono avere gli oggetti identificati
maxDim1=1500; % dimensione massima in pixel che devono avere gli oggetti identificati
maxDim2=300000;
fLayer=7; % spessore in pixel del layer che definisce l'intorno del cluster nel quale si misura il fondo da sottrarre all'intensit� del cluster stesso
sat=255; % massimo valore che pu� assumere l'intensit� del pixel (255 o 65535)
aRic=22;  % definisce in pixel met� del lato dell'area di ricerca nella quale cercare lo spostamento del centroide

dataCarlo= bfopen       % apre il filmato (hai bisogno dell'apposita libreria)

seriesC = dataCarlo{1, 1};

clear dataCarlo  % elimina la variabile dataCarlo che non serve pi� per liberare memoria

series_plane1 = seriesC{1, 1};                     % estraggo i 3 canali del primo frame dal file .lif e salvo in .tiff tutti gli stack dei 3 canali
bfsave(series_plane1,'cfp1.tiff');
series_plane2 = seriesC{2, 1};             
bfsave(series_plane2,'yfp1.tiff');                    
series_plane3 = seriesC{3, 1};
bfsave(series_plane3,'red1.tiff');
series_plane1 = seriesC{4, 1};                     
bfsave(series_plane1,'cfp2.tiff');
series_plane2 = seriesC{5, 1};             
bfsave(series_plane2,'yfp2.tiff');                      
series_plane3 = seriesC{6, 1};
bfsave(series_plane3,'red2.tiff');
series_plane1 = seriesC{7, 1};                     
bfsave(series_plane1,'cfp3.tiff');
series_plane2 = seriesC{8, 1};             
bfsave(series_plane2,'yfp3.tiff');                   
series_plane3 = seriesC{9, 1};
bfsave(series_plane3,'red3.tiff');
series_plane1 = seriesC{10, 1};                    
bfsave(series_plane1,'cfp4.tiff');
series_plane2 = seriesC{11, 1};             
bfsave(series_plane2,'yfp4.tiff');                      
series_plane3 = seriesC{12, 1};
bfsave(series_plane3,'red4.tiff');
if st > 4
    series_plane1 = seriesC{13, 1};                     
    bfsave(series_plane1,'cfp5.tiff');
    series_plane2 = seriesC{14, 1};
    bfsave(series_plane2,'yfp5.tiff');                    
    series_plane3 = seriesC{15, 1};
    bfsave(series_plane3,'red5.tiff');
end

Iy1=imread('yfp1.tiff');                           % apro le immagini dei 3 canali (come matrici)
Ic1=imread('cfp1.tiff');  
Ir1=imread('red1.tiff');
Iy2=imread('yfp2.tiff');                           
Ic2=imread('cfp2.tiff');  
Ir2=imread('red2.tiff');
Iy3=imread('yfp3.tiff');                          
Ic3=imread('cfp3.tiff');  
Ir3=imread('red3.tiff');
Iy4=imread('yfp4.tiff');                           
Ic4=imread('cfp4.tiff');  
Ir4=imread('red4.tiff');
if st > 4
    Iy5=imread('yfp5.tiff');                           
    Ic5=imread('cfp5.tiff');
    Ir5=imread('red5.tiff');
end

if st > 4
    vettImmY=cat(3,Iy1,Iy2,Iy3,Iy4,Iy5);     % unisce gli stacks in un unico array
    vettImmC=cat(3,Ic1,Ic2,Ic3,Ic4,Ic5);
    vettImmR=cat(3,Ir1,Ir2,Ir3,Ir4,Ir5);
else
    vettImmY=cat(3,Iy1,Iy2,Iy3,Iy4);
    vettImmC=cat(3,Ic1,Ic2,Ic3,Ic4);
    vettImmR=cat(3,Ir1,Ir2,Ir3,Ir4);
end

mipY = max(vettImmY, [], 3);         % fa la maximum projection degli stacks (tutti, per escludere degli stack si devono modificare i "vettImm"
mipC = max(vettImmC, [], 3);
mipR = max(vettImmR, [], 3);

% for kkk=1:512     % cicli che settano a zero i pixel con valori che sembrano corrispondere a picchi di rumore, si possono includere o meno, per verificare la presenza di questo rumore usa la funzione "imhist" sulle variabili "mip"
% for kk=1:512
% if mipY(kkk,kk) > 126 && mipY(kkk,kk) < 128
% mipY(kkk,kk)=0;
% end
% if mipR(kkk,kk) > 126 && mipR(kkk,kk) < 128
% mipR(kkk,kk)=0;
% end
% if mipC(kkk,kk) > 126 && mipC(kkk,kk) < 128
% mipC(kkk,kk)=0;
% end
% if mipC(kkk,kk) > 47 && mipC(kkk,kk) < 49
% mipC(kkk,kk)=0;
% end
% if mipY(kkk,kk) > 47 && mipY(kkk,kk) < 49
% mipY(kkk,kk)=0;
% end
% if mipR(kkk,kk) > 47 && mipR(kkk,kk) < 49
% mipR(kkk,kk)=0;
% end
% if mipC(kkk,kk) > 22 && mipC(kkk,kk) < 24
% mipC(kkk,kk)=0;
% end
% if mipR(kkk,kk) > 22 && mipR(kkk,kk) < 24
% mipR(kkk,kk)=0;
% end
% if mipY(kkk,kk) > 22 && mipY(kkk,kk) < 24
% mipY(kkk,kk)=0;
% end
% % if mipR(kkk,kk) > 32380 && mipR(kkk,kk) < 32900
% %     mipR(kkk,kk)=0;
% % end
% end
% end

% for i=1:512       % prepara l'immagine da segmentare:
% for j=1:512
% if mipR(i,j)>minR && mipC(i,j)>minC  % costruisce una nuova immagine che ha i valori di mipC solo per i pixel che hanno intensit� maggiore di minC e minR sia nel canale del ciano che del magenta, serve per eliminare parte del rumore
% mixxx(i,j)=mipC(i,j); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% else
% mixxx(i,j)=0;  % i pixel con valori minori di "minR" e/o "minC" vengono settati a 0 nella nuova immagine
% end
% end
% end
mixxx=mat2gray(mipY); %%%%%%%%%%%%%%%%%%
mixxx = imadjust(mixxx,[0;contr],[0;1]);
h = fspecial('gaussian',gaussD,gaussS);     % SMOOTHING GAUSSIANO
mixxx = filter2(h,mixxx);
mixxxbw=im2bw(mixxx,thre1);
mixxxbw=medfilt2(mixxxbw);
mixxxbw=medfilt2(mixxxbw);

mipRbw=imclose(mixxxbw, strel('disk',1));  % fa un closing dell'immagine, trasformazione che tende a uniformare i bordi degli oggetti aumentandone la dimensione
mipRbw=imfill(mipRbw,'holes');    % riempie i buchi all'interno degli oggetti
mipRbw=imdilate(mipRbw,strel('disk',1));  % dilatazioni per aumentare la dimensione degli oggetti (valuta se sono necessarie, sono legate alla soglia del thresholding)
% mipRbw=imdilate(mipRbw,strel('disk',1));
D = -bwdist(~mipRbw);     % sequenza che implementa la watershed transorm, che va effettivamente a segmentare l'immagine
D(~mipRbw) = -Inf;
LLL = watershed(D);     % controlla il risultato della segmentazione e in base a questo regola i parametri iniziali
mLLL=max(LLL(:));                   % ciclo che serve solo a settare a 0 il fondo
[a, b] = find(LLL==1);
for x = 1:length([a,b])
LLL(a(x),b(x))=0;
end
LLL=bwlabel(LLL);   % fa un relabelling dell'immagine

%%%%%%%%%%% vedi se questa sequenza � utile o meno
LLL=imdilate(LLL,strel('disk',1));   % un'altra dilatazione
mmm=im2bw(LLL,0);  % ritrasforma l'immagine in binaria e fa di nuovo il labelling
LLL=bwlabel(mmm);
%%%%%%%%%%%

for j = 1:mLLL                %ELIMINO I CLUSTER PIU' PICCOLI DI "minDim" PIXEL E PIU' GRANDI DI "maxDim"
[a, b] = find(LLL==j);
ab=[a b];
sab=size(ab);
if sab(1) < minDim
for x = 1:numel(a)
LLL(a(x),b(x))=0;
end
end
if sab(1) > maxDim1
for x = 1:numel(a)
LLL(a(x),b(x))=0;
end
end
end
LLL=imfill(LLL,'holes'); % riempie di nuovo eventuali buchi
LLL=bwlabel(LLL); % fa di nuovo il labellinge dell'immagine (necessario dopo aver eliminato degli oggetti)

uno=LLL;

mixxx=mat2gray(mipY); %%%%%%%%%%%%%%%%%%
mixxx = imadjust(mixxx,[0;contr],[0;1]);
h = fspecial('gaussian',gaussD,gaussS);     % SMOOTHING GAUSSIANO
mixxx = filter2(h,mixxx);
mixxxbw=im2bw(mixxx,thre2);
mixxxbw=medfilt2(mixxxbw);
mixxxbw=medfilt2(mixxxbw);

mipRbw=imclose(mixxxbw, strel('disk',1));
mipRbw=imfill(mipRbw,'holes');
mipRbw=imdilate(mipRbw,strel('disk',1));
% mipRbw=imdilate(mipRbw,strel('disk',1)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = -bwdist(~mipRbw);     
D(~mipRbw) = -Inf;
LLL = watershed(D);
mLLL=max(LLL(:));                   % SETTO A ZERO IL FONDO
[a, b] = find(LLL==1);
for x = 1:length([a,b])
LLL(a(x),b(x))=0;
end
LLL=bwlabel(LLL);

%%%%%%%%%%%%%%
LLL=imdilate(LLL,strel('disk',1));
mmm=im2bw(LLL,0);
LLL=bwlabel(mmm);
%%%%%%%%%%%%%%

for j = 1:mLLL                %ELIMINO I CLUSTER PIU' PICCOLI DI "minDim" PIXEL E PIU' GRANDI DI "maxDim"
[a, b] = find(LLL==j);
ab=[a b];
sab=size(ab);
if sab(1) < minDim
for x = 1:numel(a)
LLL(a(x),b(x))=0;
end
end
if sab(1) > maxDim2
for x = 1:numel(a)
LLL(a(x),b(x))=0;
end
end
end
LLL=imfill(LLL,'holes');
LLL=bwlabel(LLL);

due=LLL;

for i=1:size(uno,1)
for j=1:size(uno,2)
if uno(i,j)>0 || due(i,j)>0
tre(i,j)=1;
else
tre(i,j)=0;
end
end
end

LLL=bwlabel(tre);

vettoreFrame=cat(3,LLL);       % array che immagazzina tutti i frame segmentati, pu� essere utile per capire dopo l'analisi la qualit� della segmentazione, appesantisce un po' i calcoli
% vettFrR=cat(3,vettFrR,mipR);
% vettFrY=cat(3,vettFrY,mipY);
% vettFrC=cat(3,vettFrC,mipC);

fretMedia(1)=(mean(mean(mipR))-0.22*mean(mean(mipC)))/mean(mean(mipC));  % calcola l'attivazione su tutto il frame
RMedia(1)=mean(mean(mipR));
YMedia(1)=mean(mean(mipY));
CMedia(1)=mean(mean(mipC));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stat1=0;       
stat1 = regionprops(LLL,'centroid');           % individua i centroidi degli oggetti dell'immagine segmentata LL

for x = 1: numel(stat1)              % calcola la distanza di ogni centroide dal punto corrispondente al nucleo
    aaax=stat1(x).Centroid(1);
    bbby=stat1(x).Centroid(2);
    
    distance(1,x)=sqrt(((aaax-nXpos)^2)+((bbby-nYpos)^2));
end
    

% CICLO PER IL CALCOLO DEI VALORI DI FLUORESCENZA DEL PRIMO FRAME

   mLL=max(LLL(:));     % calcola il valore massimo assunto dalla matrice "LL"
   j=0;     
   fluoy=0;
   fluoc=0;
   fluor=0;

   for j = 0:mLL                 % loop su tutti i cluster ("0" E' IL FONDO)
       sy=0;
       sc=0;
       sr=0;
       [r, l] = find(LLL==j);
       rl=[r l];
       srl=size(rl);
       
       yma=max(rl(:,1));  % individua l'estensione massima del cluster (rettangolo nel quale � inscritto)
       ymi=min(rl(:,1));
       xma=max(rl(:,2));
       xmi=min(rl(:,2));
       
       if j ~= 0  % non � detto che questo modo di togliere il fondo sia il pi� efficace, bisognerebbe valutare se togliendo il fondo globale dell'immagine (ovvere il valore calcolato sul complementare della maschera degli oggetti segmentati) si ottengono risultati migliori
       cf=0;
       ffoy=0;
       ffoc=0;
       ffore=0;
       for ix=(xmi-fLayer):(xma+fLayer)  % definisce un intorno del cluster nel quale andare a misurare il valore del fondo da togliere all'intensit� del cluster stesso
           for iy=(ymi-fLayer):(yma+fLayer)
               if ix>0 && iy>0 && ix<513 && iy<513
               if LLL(iy,ix) ~= j
                   cf=cf+1;
                   ffoy(cf)=mipY(iy,ix);
                   ffoc(cf)=mipC(iy,ix);
                   ffore(cf)=mipR(iy,ix);
               end
               end
           end
       end
       foy(j)=sum(ffoy);
       foc(j)=sum(ffoc);
       fore(j)=sum(ffore);
       foyM(j)=foy(j)/cf; % valori misurati del fondo (per pixel) nei 3 canali
       focM(j)=foc(j)/cf;
       foreM(j)=fore(j)/cf;
       end

       i=1;
       for i = 1:srl(1)  % crea il vettore che immagazzina l'intensit� di ogni singolo pixel appartenente al cluster
           cord=rl(i,(1:2));
           sy(i)=mipY(cord(1),cord(2));
           sc(i)=mipC(cord(1),cord(2));
           sr(i)=mipR(cord(1),cord(2));
       end
       
       cy=0;  %% Non considero il contributo dei pixel saturi (metto a 0 il loro valore e non li considero quando medio)
       cc=0;
       cr=0;
       for k=1:size(sy,2)
           if sy(k)==sat  
               sy(k)=0;
               cy(k)=1;
           end
           if sc(k)==sat
               sc(k)=0;
               cc(k)=1;
           end
           if sr(k)==sat
               sr(k)=0;
               cr(k)=1;
           end
       end
       
       if j==0
           fondy=sum(sy);%/srl(1);
           fondc=sum(sc);%/srl(1);
           fondr=sum(sr);%/srl(1);
           fondsize=srl(1);
       end

       if j~=0
           fluoy(j)=sum(sy)/(srl(1)-sum(cy));  %% tolgo alla media i pixel saturi
           fluoc(j)=sum(sc)/(srl(1)-sum(cc));
           fluor(j)=sum(sr)/(srl(1)-sum(cr));
           sirl(j)=srl(1);            % immagazzina la dimensione dei cluster
       end

   end
   
   maxy=1;%max(fluoy);      % immagazzino il massimo valore di fluorescenza sui 4 canali per il primo frame
   maxc=1;%max(fluoc);    % user� questi valori per la normalizzazione in tutti i frame
   maxr=1;%max(fluor);
   
   fondoy(1)=fondy/maxy;
   fondoc(1)=fondc/maxc;
   fondor(1)=fondr/maxr;
   mfondoy(1)=fondoy(1)/fondsize;           % tengo traccia del fondo medio per pixel
   mfondoc(1)=fondoc(1)/fondsize;
   mfondor(1)=fondor(1)/fondsize;

   for jj = 1:mLL   
       ffluoy(1,jj)=fluoy(jj);%-foyM(jj);  %AGGIORNA LE MATRICI DELLE FLUORESCENZE NORMALIZZATE sottraendo l'autofluorescenza
       ffluoc(1,jj)=fluoc(jj);%-focM(jj);                      
       ffluor(1,jj)=fluor(jj);%-foreM(jj);     
       ssirl(1,jj)=sirl(jj);
   end
   
  
ftoty(1)=sum(sum(mipY))/maxy;%(262144*maxy);         % SALVA L'INTENSITA' TOTALE DEL FRAME (oppure intensit� media per pixel)
ftotc(1)=sum(sum(mipC))/maxc;%(262144*maxc);
ftotch(1)=sum(sum(mipR))/maxr;%(262144*maxr);

fcly(1)=(ftoty(1)-fondoy(1));    % intensit� globale dei pixel corrispondenti agli oggetti segmentati
fclc(1)=(ftotc(1)-fondoc(1));
fclr(1)=(ftotr(1)-fondor(1));

% for cn1=r1x1:r1x2   % calcolo le intensita` medie in zone del campo non coperte da vescicole
%     for cn2=r1y1:r1y2
%         ret1R(cn2-r1y1+1,cn1-r1x1+1)=mipR(cn2,cn1);
%         ret1C(cn2-r1y1+1,cn1-r1x1+1)=mipC(cn2,cn1);
%         ret1Y(cn2-r1y1+1,cn1-r1x1+1)=mipY(cn2,cn1);
%     end
% end
% clear cn1 cn2
% for cn1=r2x1:r2x2
%     for cn2=r2y1:r2y2
%         ret2R(cn2-r2y1+1,cn1-r2x1+1)=mipR(cn2,cn1);
%         ret2C(cn2-r2y1+1,cn1-r2x1+1)=mipC(cn2,cn1);
%         ret2Y(cn2-r2y1+1,cn1-r2x1+1)=mipY(cn2,cn1);
%     end
% end
% clear cn1 cn2
% for cn1=r3x1:r3x2
%     for cn2=r3y1:r3y2
%         ret3R(cn2-r3y1+1,cn1-r3x1+1)=mipR(cn2,cn1);
%         ret3C(cn2-r3y1+1,cn1-r3x1+1)=mipC(cn2,cn1);
%         ret3Y(cn2-r3y1+1,cn1-r3x1+1)=mipY(cn2,cn1);
%     end
% end
% clear cn1 cn2
% for cn1=r4x1:r4x2
%     for cn2=r4y1:r4y2
%         ret4R(cn2-r4y1+1,cn1-r4x1+1)=mipR(cn2,cn1);
%         ret4C(cn2-r4y1+1,cn1-r4x1+1)=mipC(cn2,cn1);
%         ret4Y(cn2-r4y1+1,cn1-r4x1+1)=mipY(cn2,cn1);
%     end
% end
% RET1R(1)=mean(mean(ret1R));
% RET1C(1)=mean(mean(ret1C));
% RET1Y(1)=mean(mean(ret1Y));
% RET2R(1)=mean(mean(ret2R));
% RET2C(1)=mean(mean(ret2C));
% RET2Y(1)=mean(mean(ret2Y));
% RET3R(1)=mean(mean(ret3R));
% RET3C(1)=mean(mean(ret3C));
% RET3Y(1)=mean(mean(ret3Y));
% RET4R(1)=mean(mean(ret4R));
% RET4C(1)=mean(mean(ret4C));
% RET4Y(1)=mean(mean(ret4Y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


delete yfp1.tiff                       % elimino le immagini create precedentemente in modo da poterle salvare nuovamente
delete cfp1.tiff
delete red1.tiff
delete yfp2.tiff                       % elimino le immagini create precedentemente in modo da poterle salvare nuovamente
delete cfp2.tiff
delete red2.tiff
delete yfp3.tiff                       % elimino le immagini create precedentemente in modo da poterle salvare nuovamente
delete cfp3.tiff
delete red3.tiff
delete yfp4.tiff                       % elimino le immagini create precedentemente in modo da poterle salvare nuovamente
delete cfp4.tiff
delete red4.tiff
if st>4
    delete yfp5.tiff                       % elimino le immagini create precedentemente in modo da poterle salvare nuovamente
    delete cfp5.tiff
    delete red5.tiff
end

LL=LLL;

mLL=max(LL(:));

for K=2:size(seriesC,1)/(3*st) % inizia il loop su tutti i frame del filmato

%      if K>18 %&& K<28
%           st=5; % numero stacks (4 o 5 nel caso siano di pi� o di meno ci sono da fare piccole modifiche)
%  minC=35; 
%  minR=35; % sono i valori minimi di intensit� del pixel mecessari per considerare segnale reale il dato (vedi punto di utilizzo per maggiori info)
%  gaussD=[3 3]; % dimensione del filtro gaussiano
%  gaussS=0.95; % deviazione standard del filtro gaussiano
%  thre1=0.28; % soglia del thresholding che trasforma l'immagine da segmentare in binaria
%  thre2=0.8;
%  minDim=60; % dimensione minima in pixel che devono avere gli oggetti identificati
%  maxDim1=2000; % dimensione massima in pixel che devono avere gli oggetti identificati
%  maxDim2=300000;
%  fLayer=7; % spessore in pixel del layer che definisce l'intorno del cluster nel quale si misura il fondo da sottrarre all'intensit� del cluster stesso
%  sat=255; % massimo valore che pu� assumere l'intensit� del pixel (255 o 65535)
%  aRic=22; % definisce in pixel met� del lato dell'area di ricerca nella quale cercare lo spostamento del centroide
%       end
%     if K>28 && K<54
%         st=5; % numero stacks (4 o 5 nel caso siano di pi� o di meno ci sono da fare piccole modifiche)
%         minC=35;
%         minR=35; % sono i valori minimi di intensit� del pixel mecessari per considerare segnale reale il dato (vedi punto di utilizzo per maggiori info)
%         gaussD=[3 3]; % dimensione del filtro gaussiano
%         gaussS=0.95; % deviazione standard del filtro gaussiano
%         thre=0.32; % soglia del thresholding che trasforma l'immagine da segmentare in binaria
%         minDim=40; % dimensione minima in pixel che devono avere gli oggetti identificati
%         maxDim=50000; % dimensione massima in pixel che devono avere gli oggetti identificati
%         fLayer=7; % spessore in pixel del layer che definisce l'intorno del cluster nel quale si misura il fondo da sottrarre all'intensit� del cluster stesso
%         sat=255; % massimo valore che pu� assumere l'intensit� del pixel (255 o 65535)
%         aRic=22; % definisce in pixel met� del lato dell'area di ricerca nella quale cercare lo spostamento del centroide
%     end
%     if K<54
%         st=5; % numero stacks (4 o 5 nel caso siano di pi� o di meno ci sono da fare piccole modifiche)
%         minC=35;
%         minR=35; % sono i valori minimi di intensit� del pixel mecessari per considerare segnale reale il dato (vedi punto di utilizzo per maggiori info)
%         gaussD=[3 3]; % dimensione del filtro gaussiano
%         gaussS=0.95; % deviazione standard del filtro gaussiano
%         thre=0.25; % soglia del thresholding che trasforma l'immagine da segmentare in binaria
%         minDim=40; % dimensione minima in pixel che devono avere gli oggetti identificati
%         maxDim=5000; % dimensione massima in pixel che devono avere gli oggetti identificati
%         fLayer=7; % spessore in pixel del layer che definisce l'intorno del cluster nel quale si misura il fondo da sottrarre all'intensit� del cluster stesso
%         sat=255; % massimo valore che pu� assumere l'intensit� del pixel (255 o 65535)
%         aRic=22; % definisce in pixel met� del lato dell'area di ricerca nella quale cercare lo spostamento del centroide
%     end

    
    
    
clear counter
counter=ones(1,mLL);  % definisco il vettore che terr� conto delle corrispondeze tra i cluster di un frame e quelli del successivo

series_plane1 = seriesC{(K-1)*(3*st)+1, 1};                     % estraggo i 3 canali del frame dal file .lif
bfsave(series_plane1,'cfp1.tiff');
series_plane2 = seriesC{(K-1)*(3*st)+2, 1};             
bfsave(series_plane2,'yfp1.tiff');                    % K � l'indice del frame 
series_plane3 = seriesC{(K-1)*(3*st)+3, 1};
bfsave(series_plane3,'red1.tiff');
series_plane1 = seriesC{(K-1)*(3*st)+4, 1};                     
bfsave(series_plane1,'cfp2.tiff');
series_plane2 = seriesC{(K-1)*(3*st)+5, 1};             
bfsave(series_plane2,'yfp2.tiff');                   
series_plane3 = seriesC{(K-1)*(3*st)+6, 1};
bfsave(series_plane3,'red2.tiff');
series_plane1 = seriesC{(K-1)*(3*st)+7, 1};                     
bfsave(series_plane1,'cfp3.tiff');
series_plane2 = seriesC{(K-1)*(3*st)+8, 1};             
bfsave(series_plane2,'yfp3.tiff');                   
series_plane3 = seriesC{(K-1)*(3*st)+9, 1};
bfsave(series_plane3,'red3.tiff');
series_plane1 = seriesC{(K-1)*(3*st)+10, 1};                     
bfsave(series_plane1,'cfp4.tiff');
series_plane2 = seriesC{(K-1)*(3*st)+11, 1};             
bfsave(series_plane2,'yfp4.tiff');                     
series_plane3 = seriesC{(K-1)*(3*st)+12, 1};
bfsave(series_plane3,'red4.tiff');
if st > 4
    series_plane1 = seriesC{(K-1)*(3*st)+13, 1};                     
    bfsave(series_plane1,'cfp5.tiff');
    series_plane2 = seriesC{(K-1)*(3*st)+14, 1};
    bfsave(series_plane2,'yfp5.tiff');                    
    series_plane3 = seriesC{(K-1)*(3*st)+15, 1};
    bfsave(series_plane3,'red5.tiff');
end

Iy1=imread('yfp1.tiff');                           % apro le immagini dei 3 canali (come matrici)
Ic1=imread('cfp1.tiff');  
Ir1=imread('red1.tiff');
Iy2=imread('yfp2.tiff');                           
Ic2=imread('cfp2.tiff');  
Ir2=imread('red2.tiff');
Iy3=imread('yfp3.tiff');                           
Ic3=imread('cfp3.tiff');  
Ir3=imread('red3.tiff');
Iy4=imread('yfp4.tiff');                           
Ic4=imread('cfp4.tiff');  
Ir4=imread('red4.tiff');
if st > 4
    Iy5=imread('yfp5.tiff');                           
    Ic5=imread('cfp5.tiff');
    Ir5=imread('red5.tiff');
end

clear vettImmY vettImmC vettImmR

if st > 4
    vettImmY=cat(3,Iy1,Iy2,Iy3,Iy4,Iy5);
    vettImmC=cat(3,Ic1,Ic2,Ic3,Ic4,Ic5);
    vettImmR=cat(3,Ir1,Ir2,Ir3,Ir4,Ir5);
else
    vettImmY=cat(3,Iy1,Iy2,Iy3,Iy4);
    vettImmC=cat(3,Ic1,Ic2,Ic3,Ic4);
    vettImmR=cat(3,Ir1,Ir2,Ir3,Ir4);
end

mipY = max(vettImmY, [], 3);
mipC = max(vettImmC, [], 3);
mipR = max(vettImmR, [], 3);

% for kkk=1:512
% for kk=1:512
% if mipY(kkk,kk) > 126 && mipY(kkk,kk) < 128
% mipY(kkk,kk)=0;
% end
% if mipR(kkk,kk) > 126 && mipR(kkk,kk) < 128
% mipR(kkk,kk)=0;
% end
% if mipC(kkk,kk) > 126 && mipC(kkk,kk) < 128
% mipC(kkk,kk)=0;
% end
% if mipC(kkk,kk) > 47 && mipC(kkk,kk) < 49
% mipC(kkk,kk)=0;
% end
% if mipY(kkk,kk) > 47 && mipY(kkk,kk) < 49
% mipY(kkk,kk)=0;
% end
% if mipR(kkk,kk) > 47 && mipR(kkk,kk) < 49
% mipR(kkk,kk)=0;
% end
% if mipC(kkk,kk) > 22 && mipC(kkk,kk) < 24
% mipC(kkk,kk)=0;
% end
% if mipR(kkk,kk) > 22 && mipR(kkk,kk) < 24
% mipR(kkk,kk)=0;
% end
% if mipY(kkk,kk) > 22 && mipY(kkk,kk) < 24
% mipY(kkk,kk)=0;
% end
% % if mipR(kkk,kk) > 32380 && mipR(kkk,kk) < 32900
% %     mipR(kkk,kk)=0;
% % end
% end
% end

% for i=1:512
% for j=1:512
% if mipR(i,j)>minR && mipC(i,j)>minC
% mixxx(i,j)=mipR(i,j); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% else
% mixxx(i,j)=0;
% end
% end
% end
mixxx=mat2gray(mipY); %%%%%%%%%%%%%%%%%%
mixxx = imadjust(mixxx,[0;contr],[0;1]);
h = fspecial('gaussian',gaussD,gaussS);     % SMOOTHING GAUSSIANO
mixxx = filter2(h,mixxx);
mixxxbw=im2bw(mixxx,thre1);
mixxxbw=medfilt2(mixxxbw);
mixxxbw=medfilt2(mixxxbw);

mipRbw=imclose(mixxxbw, strel('disk',1));
mipRbw=imfill(mipRbw,'holes');
mipRbw=imdilate(mipRbw,strel('disk',1));
% mipRbw=imdilate(mipRbw,strel('disk',1)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = -bwdist(~mipRbw);     
D(~mipRbw) = -Inf;
LLL = watershed(D);
mLLL=max(LLL(:));                   % SETTO A ZERO IL FONDO
[a, b] = find(LLL==1);
for x = 1:length([a,b])
LLL(a(x),b(x))=0;
end
LLL=bwlabel(LLL);

%%%%%%%%%%%%%%
LLL=imdilate(LLL,strel('disk',1));
mmm=im2bw(LLL,0);
LLL=bwlabel(mmm);
%%%%%%%%%%%%%%

for j = 1:mLLL                %ELIMINO I CLUSTER PIU' PICCOLI DI "minDim" PIXEL E PIU' GRANDI DI "maxDim"
[a, b] = find(LLL==j);
ab=[a b];
sab=size(ab);
if sab(1) < minDim
for x = 1:numel(a)
LLL(a(x),b(x))=0;
end
end
if sab(1) > maxDim1
for x = 1:numel(a)
LLL(a(x),b(x))=0;
end
end
end
LLL=imfill(LLL,'holes');
LLL=bwlabel(LLL);

uno=LLL;

mixxx=mat2gray(mipY); %%%%%%%%%%%%%%%%%%
mixxx = imadjust(mixxx,[0;contr],[0;1]);
h = fspecial('gaussian',gaussD,gaussS);     % SMOOTHING GAUSSIANO
mixxx = filter2(h,mixxx);
mixxxbw=im2bw(mixxx,thre2);
mixxxbw=medfilt2(mixxxbw);
mixxxbw=medfilt2(mixxxbw);

mipRbw=imclose(mixxxbw, strel('disk',1));
mipRbw=imfill(mipRbw,'holes');
mipRbw=imdilate(mipRbw,strel('disk',1));
% mipRbw=imdilate(mipRbw,strel('disk',1)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = -bwdist(~mipRbw);     
D(~mipRbw) = -Inf;
LLL = watershed(D);
mLLL=max(LLL(:));                   % SETTO A ZERO IL FONDO
[a, b] = find(LLL==1);
for x = 1:length([a,b])
LLL(a(x),b(x))=0;
end
LLL=bwlabel(LLL);

%%%%%%%%%%%%%%
LLL=imdilate(LLL,strel('disk',1));
mmm=im2bw(LLL,0);
LLL=bwlabel(mmm);
%%%%%%%%%%%%%%

for j = 1:mLLL                %ELIMINO I CLUSTER PIU' PICCOLI DI "minDim" PIXEL E PIU' GRANDI DI "maxDim"
[a, b] = find(LLL==j);
ab=[a b];
sab=size(ab);
if sab(1) < minDim
for x = 1:numel(a)
LLL(a(x),b(x))=0;
end
end
if sab(1) > maxDim2
for x = 1:numel(a)
LLL(a(x),b(x))=0;
end
end
end
LLL=imfill(LLL,'holes');
LLL=bwlabel(LLL);

due=LLL;

for i=1:size(uno,1)
for j=1:size(uno,2)
if uno(i,j)>0 || due(i,j)>0
tre(i,j)=1;
else
tre(i,j)=0;
end
end
end

LLL=bwlabel(tre);

LL=LLL;

mLL=max(LL(:));

vettoreFrame=cat(3,vettoreFrame,LL);    % aggiorna il vettore che immagazzina i frame
% vettFrR=cat(3,vettFrR,mipR);
% vettFrY=cat(3,vettFrY,mipY);
% vettFrC=cat(3,vettFrC,mipC);

fretMedia(K)=(mean(mean(mipR))-0.22*mean(mean(mipC)))/mean(mean(mipC));
RMedia(K)=mean(mean(mipR));
YMedia(K)=mean(mean(mipY));
CMedia(K)=mean(mean(mipC));

stat1=0;       
stat1 = regionprops(LL,'centroid');           % individua i centroidi degli oggetti dell'immagine segmentata LL

for x = 1: numel(stat1)              % calcola la distanza di ogni centroide dal punto corrispondente al nucleo
    aaax=stat1(x).Centroid(1);
    bbby=stat1(x).Centroid(2);
    
    distance(K,x)=sqrt(((aaax-nXpos)^2)+((bbby-nYpos)^2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CICLO PER IL CALCOLO DEI VALORI DI FLUORESCENZA DEI SINGOLI FRAME

   j=0;     
   fluoy=0;
   fluoc=0;
   fluor=0;

   for j = 0:mLL                 % loop su tutti i cluster ("0" E' IL FONDO)
       sy=0;
       sc=0;
       sr=0;
       [r, l] = find(LL==j);
       rl=[r l];
       srl=size(rl);
       
       yma=max(rl(:,1));
       ymi=min(rl(:,1));
       xma=max(rl(:,2));
       xmi=min(rl(:,2));
       
       if j ~= 0
       cf=0;
       ffoy=0;
       ffoc=0;
       ffore=0;
       for ix=(xmi-fLayer):(xma+fLayer)
           for iy=(ymi-fLayer):(yma+fLayer)
               if ix>0 && iy>0 && ix<513 && iy<513
               if LLL(iy,ix) ~= j
                   cf=cf+1;
                   ffoy(cf)=mipY(iy,ix);
                   ffoc(cf)=mipC(iy,ix);
                   ffore(cf)=mipR(iy,ix);
               end
               end
           end
       end
       foy(j)=sum(ffoy);
       foc(j)=sum(ffoc);
       fore(j)=sum(ffore);
       foyM(j)=foy(j)/cf;
       focM(j)=foc(j)/cf;
       foreM(j)=fore(j)/cf;
       end

       i=1;
       for i = 1:srl(1)
           cord=rl(i,(1:2));
           sy(i)=mipY(cord(1),cord(2));
           sc(i)=mipC(cord(1),cord(2));
           sr(i)=mipR(cord(1),cord(2));
       end
       
       cy=0;  %% Non considero il contributo dei pixel saturi (metto a 0 il loro valore e non li considero quando medio)
       cc=0;
       cr=0;
       for k=1:size(sy,2)
           if sy(k)==sat  %% Settare il valore giusto per la saturazione !!!!
               sy(k)=0;
               cy(k)=1;
           end
           if sc(k)==sat
               sc(k)=0;
               cc(k)=1;
           end
           if sr(k)==sat
               sr(k)=0;
               cr(k)=1;
           end
       end
       
       if j==0
           fondy=sum(sy);%/srl(1);
           fondc=sum(sc);%/srl(1);
           fondr=sum(sr);%/srl(1);
           fondsize=srl(1);
       end

       if j~=0
           fluoy(j)=sum(sy)/(srl(1)-sum(cy));  %% tolgo alla media i pixel saturi
           fluoc(j)=sum(sc)/(srl(1)-sum(cc));
           fluor(j)=sum(sr)/(srl(1)-sum(cr));
           sirl(j)=srl(1);            % immagazzina la dimensione dei cluster
       end

   end
   
   maxy=1;%max(fluoy);      % immagazzino il massimo valore di fluorescenza sui 4 canali per il primo frame
   maxc=1;%max(fluoc);    % user� questi valori per la normalizzazione in tutti i frame
   maxr=1;%max(fluor);
   
   fondoy(K)=fondy/maxy;
   fondoc(K)=fondc/maxc;
   fondor(K)=fondr/maxr;
   mfondoy(K)=fondoy(K)/fondsize;           % tengo traccia del fondo medio per pixel
   mfondoc(K)=fondoc(K)/fondsize;
   mfondor(K)=fondor(K)/fondsize;

   for jj = 1:mLL   
       ffluoy(K,jj)=fluoy(jj);%-foyM(jj);%/max(fluoy)-mfondoy(k);%sirl(jj)*mfondoy(1);   %AGGIORNA LE MATRICI DELLE FLUORESCENZE NORMALIZZATE sottraendo l'autofluorescenza
       ffluoc(K,jj)=fluoc(jj);%-focM(jj);%/max(fluoc)-mfondoc(k);%sirl(jj)*mfondoc(1);                       
       ffluor(K,jj)=fluor(jj);%-foreM(jj);%/max(fluor)-mfondor(k);%sirl(jj)*mfondor(1);      
       ssirl(K,jj)=sirl(jj);
   end   
  
ftoty(K)=sum(sum(mipY))/maxy;%(262144*maxy);         % SALVA L'INTENSITA' TOTALE DEL FRAME (intensit� media per pixel)
ftotc(K)=sum(sum(mipC))/maxc;%(262144*maxc);
ftotr(K)=sum(sum(mipR))/maxr;%(262144*maxr);

fcly(K)=(ftoty(K)-fondoy(K));
fclc(K)=(ftotc(K)-fondoc(K));
fclr(K)=(ftotr(K)-fondor(K));

% for cn1=r1x1:r1x2  % calcolo le intensita` medie in zone del campo non coperte da vescicole
%     for cn2=r1y1:r1y2
%         ret1R(cn2-r1y1+1,cn1-r1x1+1)=mipR(cn2,cn1);
%         ret1C(cn2-r1y1+1,cn1-r1x1+1)=mipC(cn2,cn1);
%         ret1Y(cn2-r1y1+1,cn1-r1x1+1)=mipY(cn2,cn1);
%     end
% end
% clear cn1 cn2
% for cn1=r2x1:r2x2
%     for cn2=r2y1:r2y2
%         ret2R(cn2-r2y1+1,cn1-r2x1+1)=mipR(cn2,cn1);
%         ret2C(cn2-r2y1+1,cn1-r2x1+1)=mipC(cn2,cn1);
%         ret2Y(cn2-r2y1+1,cn1-r2x1+1)=mipY(cn2,cn1);
%     end
% end
% clear cn1 cn2
% for cn1=r3x1:r3x2
%     for cn2=r3y1:r3y2
%         ret3R(cn2-r3y1+1,cn1-r3x1+1)=mipR(cn2,cn1);
%         ret3C(cn2-r3y1+1,cn1-r3x1+1)=mipC(cn2,cn1);
%         ret3Y(cn2-r3y1+1,cn1-r3x1+1)=mipY(cn2,cn1);
%     end
% end
% clear cn1 cn2
% for cn1=r4x1:r4x2
%     for cn2=r4y1:r4y2
%         ret4R(cn2-r4y1+1,cn1-r4x1+1)=mipR(cn2,cn1);
%         ret4C(cn2-r4y1+1,cn1-r4x1+1)=mipC(cn2,cn1);
%         ret4Y(cn2-r4y1+1,cn1-r4x1+1)=mipY(cn2,cn1);
%     end
% end
% RET1R(K)=mean(mean(ret1R));
% RET1C(K)=mean(mean(ret1C));
% RET1Y(K)=mean(mean(ret1Y));
% RET2R(K)=mean(mean(ret2R));
% RET2C(K)=mean(mean(ret2C));
% RET2Y(K)=mean(mean(ret2Y));
% RET3R(K)=mean(mean(ret3R));
% RET3C(K)=mean(mean(ret3C));
% RET3Y(K)=mean(mean(ret3Y));
% RET4R(K)=mean(mean(ret4R));
% RET4C(K)=mean(mean(ret4C));
% RET4Y(K)=mean(mean(ret4Y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


delete yfp1.tiff                       % elimino le immagini create precedentemente in modo da poterle salvare nuovamente
delete cfp1.tiff
delete red1.tiff
delete yfp2.tiff                       % elimino le immagini create precedentemente in modo da poterle salvare nuovamente
delete cfp2.tiff
delete red2.tiff
delete yfp3.tiff                       % elimino le immagini create precedentemente in modo da poterle salvare nuovamente
delete cfp3.tiff
delete red3.tiff
delete yfp4.tiff                       % elimino le immagini create precedentemente in modo da poterle salvare nuovamente
delete cfp4.tiff
delete red4.tiff
if st>4
    delete yfp5.tiff                       % elimino le immagini create precedentemente in modo da poterle salvare nuovamente
    delete cfp5.tiff
    delete red5.tiff
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


frameNumber=size(ffluoc,1);
tempo=(0:1:frameNumber*1);

% ANDAMENTO NUMERO DI VESCICOLE NEL TEMPO

vesNum=zeros(1,size(ssirl,1));
for g=1:size(ssirl,1)
for gg=1:size(ssirl,2)
if ssirl(g,gg)~=0
vesNum(g)=vesNum(g)+1;
end
end
end
figure, title('Number of vesicles (for square micrometer) vs time'), hold on, plot(tempo(1:size(vesNum,2))*0.5,vesNum(:)/1024,'r')

% ANDAMENTO TAGLIA VESCICOLE NEL TEMPO

for j=1:size(ssirl,2)        
for i=1:size(ssirl,1)
if ssirl(i,j)~= 0
gg(j)=i;  %% tiene conto di quando compare la vescicola j-esima
break
end
end
end
ssirl2=zeros(size(ssirl,1),size(ssirl,2));
for j=1:size(ssirl,2)
ssirl2(1:size(ssirl,1)-gg(j)+1,j)=ssirl(gg(j):size(ssirl,1),j);
end
clear Mssirl SDssirl
for i=1:size(ssirl,1)
    clear taglie
    taglie=ssirl(i,:);
    for k=size(taglie,2):-1:1
                if taglie(k)==0
                   taglie(k)=[];
                end     
    end
    Mssirl(i)=mean(taglie);
    SDssirl(i)=std(taglie);
    Totssirl(i)=sum(taglie);
end
for i=1:size(ssirl2,1)
    clear taglie
    taglie=ssirl2(i,:);
    for k=size(taglie,2):-1:1
                if taglie(k)==0
                   taglie(k)=[];
                end     
    end
    Mssirl2(i)=mean(taglie);
end
figure, title('Vesicle mean size (nanometers) vs time'), hold on, plot(tempo(1:size(vesNum,2))*0.5,sqrt(Mssirl(:)/pi)*62.89,'r')
figure, title('Vesicle mean size (nanometers) vs time SYNCH'), hold on, plot(tempo(1:size(vesNum,2))*0.5,sqrt(Mssirl2(:)/pi)*62.89,'r')
figure, title('SD on vesicle mean size (nanometers) vs time'), hold on, plot(tempo(1:size(vesNum,2))*0.5,sqrt(SDssirl(:)/pi)*62.89,'r')
figure, title('Total vesicle size (for square micrometer) vs time'), hold on, plot(tempo(1:size(vesNum,2))*0.5,Totssirl(:)/262144,'r')

% DISTRIBUZIONE DELLE TAGLIE A FISSATO FRAME

% for aaa=1:5:65
aaa=3; % numero frame
frameSizeDistr=ssirl(aaa,:);
for aa=size(frameSizeDistr,2):-1:1
    if frameSizeDistr(aa)==0
        frameSizeDistr(aa)=[];
    end
end
figure, title('Size distribution (pixel) at frame 37'), hold on, histogram(frameSizeDistr,100)
figure, title('Size distribution (pixel) at frame 37 ZOOM'), hold on, histogram(frameSizeDistr,linspace(0,400,60))
% end

% CREA VIDEO CON LE IMMAGINI SEGMENTATE

% writerObj = VideoWriter('transferrin647_20ugml2_4_SEGMENTED_ON_BLUE');
% open(writerObj);           %% crea un filmato con le immagini segmentate
% for k = 1:size(vettoreFrame,3)
% frame = mat2gray(vettoreFrame(:,:,k));
% writeVideo(writerObj,frame);
% end
% close(writerObj);

% TRAIETTORIA SINGOLA VESCICOLA

ncell=28;
% figure;
% plot(tempo(1:frameNumber),ffluor(:,ncell),'r')
% hold on
% plot(tempo(1:frameNumber),ffluoy(:,ncell),'y')
% hold on
% % figure;
% plot(tempo(1:frameNumber),ffluoc(:,ncell),'c')
% hold off
% 

% ANALISI DELL'ATTIVAZIONE FRET

s=size(ffluoc);
for i=1:s(1)
    for j=1:s(2)
        if isnan(ffluoc(i,j))
            ffluoc(i,j)=0;
        end
        if isnan(ffluor(i,j))
            ffluor(i,j)=0;
        end
        
        if ffluor(i,j) ~= 0
            Mfret(i,j)=(ffluor(i,j)-0.22*ffluoc(i,j))/(ffluoc(i,j));%*sqrt(ssirl(i,j)/pi)*62.89);
            if ffluor(i,j)/ffluoc(i,j)<2.5
                intY(i,j)=ffluor(i,j)*ssirl(i,j);
                intC(i,j)=ffluoc(i,j)*ssirl(i,j);
            else
                intY(i,j)=0;
                intC(i,j)=0;
            end
        else
            Mfret(i,j)=0;
            intY(i,j)=0;
            intC(i,j)=0;
        end
    end
end
for j=1:size(Mfret,2)
    for i=1:size(Mfret,1)
        if Mfret(i,j)~= 0
            gg(j)=i;  %% tiene conto di quando compare la vescicola j-esima
            break
        end
    end
end
Mfret2=zeros(size(Mfret,1),size(Mfret,2));
for j=1:size(Mfret,2)
    Mfret2(1:size(Mfret,1)-gg(j)+1,j)=Mfret(gg(j):size(Mfret,1),j);
end
sumY=sum(ffluor,2); % somme della fluorescenza di tutti i cluster, senza cutoff
sumC=sum(ffluoc,2);
for i=1:size(Mfret,1)
    for j=1:size(Mfret,2)
        if Mfret(i,j)<0 || Mfret(i,j)>2.5
            Mfret(i,j)=0;
        end
        if Mfret2(i,j)<0 || Mfret2(i,j)>2.5
            Mfret2(i,j)=0;
        end
    end
end
SMfret=sum(Mfret,2); % somma dell'attivazione dei cluster
for k=1:size(Mfret,1)
    ugo=Mfret(k,:);
    ugo=ugo';
    for kk=size(ugo,1):-1:1
        if ugo(kk)==0
            ugo(kk)=[];
        end
    end
    %     SSMfret(k)=sum(ugo);
    SDMfret(k)=std(ugo);  % deviazione standard
end
SMfret2=sum(Mfret2,2); % somma dell'attivazione dei cluster
for k=1:size(Mfret2,1)
    ugo=Mfret2(k,:);
    ugo=ugo';
    for kk=size(ugo,1):-1:1
        if ugo(kk)==0
            ugo(kk)=[];
        end
    end
    %     SSMfret2(k)=sum(ugo);
    SDMfret2(k)=std(ugo); % deviazione standard
end
SintY=sum(intY,2); % somme della fluorescenza assoluta di tutti i cluster
SintC=sum(intC,2);
for k=1:size(seriesC,1)/(3*st) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ratio(k)=(sumY(k)-0.22*sumC(k))/sumC(k); % attivazione globale ottenuta dalla somma delle fluorescenze medie dei cluster (meglio "Gfret")
    fret(k)=SMfret(k)/nnz(Mfret(k,:)); % attivazioni medie
    fret2(k)=SMfret2(k)/nnz(Mfret2(k,:));
    Gfret(k)=(SintY(k)-0.22*SintC(k))/SintC(k); % attivazione globale (solo pixel corrispondenti ai cluster)
end
GfretRate(1)=0;
fretRate(1)=0;
fretRate2(1)=0;
for k=2:size(seriesC,1)/(3*st)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fretRate(k)=2*(fret(k)-fret(k-1));
    fretRate2(k)=2*(fret2(k)-fret2(k-1));
    GfretRate(k)=2*(Gfret(k)-Gfret(k-1));
end

figure, title('Yfp - Cfp ratio (all vesicles)'), hold on,
plot(tempo(1:frameNumber)*0.5,ratio/ratio(1),'g')
figure, title('Mean Yfp - Cfp ratio (FRET) '), hold on,
plot(tempo(1:frameNumber)*0.5,fret/fret(1),'b'), hold on, xlabel('time (s)'), hold on, ylabel('FRET activation'), hold on,
errorbar(tempo(1:frameNumber)*0.5,fret/fret(1),SDMfret,'gx')%, axis([-2 35 0 3])
figure, title('Mean FRET rate '), hold on,
plot(tempo(1:frameNumber)*0.5,fretRate/fret(1),'g')
figure, title('Mean Yfp - Cfp ratio (FRET) Sync '), hold on,
plot(tempo(1:frameNumber)*0.5,fret2/fret2(1),'b'), axis([0 35 0 2.5]), hold on,
errorbar(tempo(1:frameNumber)*0.5,fret2/fret2(1),SDMfret2,'gx')
figure, title('Mean FRET rate Sync '), hold on,
plot(tempo(1:frameNumber)*0.5,fretRate2/fret2(1),'g')
figure, title('Global Yfp - Cfp ratio (FRET)'), hold on, xlabel('time (s)'), hold on, ylabel('FRET activation'), hold on,
plot(tempo(1:frameNumber)*0.5,Gfret/Gfret(1),'g')
figure, title('Global FRET rate '), hold on,
plot(tempo(1:frameNumber)*0.5,GfretRate/Gfret(1),'g')

figure, title('Fret LEICA '), hold on,
plot(tempo(1:frameNumber)*0.5,fretMedia,'g')

% DISTRIBUZIONE DELL'ATTIVAZIONE A FISSATO FRAME

aaa=3;
frameFretDistr=Mfret(aaa,:);
for aa=size(frameFretDistr,2):-1:1
    if frameFretDistr(aa)==0
        frameFretDistr(aa)=[];
    end
end
figure, title('Fret distribution at frame 30'), hold on, histogram(frameFretDistr,30)
% figure, title('Size distribution (pixel) at frame 37 ZOOM'), hold on, histogram(frameSizeDistr,linspace(0,400,60))


% SCATTER PLOT FRET - DIMENSIONI (TUTTE LE VESCICOLE DI TUTTI I FRAME)

clear tgl frt fy fc tgl2 tgl3
tgl=ssirl(1,:);
tgl2=ssirl(1,:);
tgl3=ssirl(1,:);
frt=Mfret(1,:);
fy=ffluoy(1,:);
fc=ffluoc(1,:);
for k=2:size(Mfret,1)
tgl=cat(2,tgl,ssirl(k,:));
tgl2=cat(2,tgl2,ssirl(k,:));
tgl3=cat(2,tgl3,ssirl(k,:));
frt=cat(2,frt,Mfret(k,:));
fy=cat(2,fy,ffluoy(k,:));
fc=cat(2,fc,ffluoc(k,:));
end
for k=size(tgl,2):-1:1
if frt(k)==0 || tgl(k)>1200
frt(k)=[];
tgl(k)=[];
end
if fy(k)<1 || tgl2(k)>1200
fy(k)=[];
tgl2(k)=[];
end
if fc(k)<1 || tgl3(k)>1200
fc(k)=[];
tgl3(k)=[];
end
end
figure, title('Vescicle radius (nm) vs FRET'), hold on, scatter(frt(:),tgl(:), 'b','.'), hold on, xlabel('FRET'), hold on, ylabel('Radius (nm)')
figure, title('Vescicle radius (nm) vs cyan'), hold on, scatter(fc(:),tgl3(:), 'b','.'), hold on, xlabel('cyan'), hold on, ylabel('Radius (nm)')
figure, title('Vescicle radius (nm) vs yfp'), hold on, scatter(fy(:),tgl2(:), 'b','.'), hold on, xlabel('yfp'), hold on, ylabel('Radius (nm)')

%%% scatter plot di attivazione fret e fluorescenza nel CFP

citoR=RET1R';
citoR=cat(2,citoR,RET2R');
citoR=cat(2,citoR,RET3R');
citoR=cat(2,citoR,RET4R');
citoC=RET1C';
citoC=cat(2,citoC,RET2C');
citoC=cat(2,citoC,RET3C');
citoC=cat(2,citoC,RET4C');
citoY=RET1Y';
citoY=cat(2,citoY,RET2Y');
citoY=cat(2,citoY,RET3Y');
citoY=cat(2,citoY,RET4Y');
clear citR citC citY
n=64;
citY=citoY(n,:);
citC=citoC(n,:);
citR=citoR(n,:);
for k=n:n+4
citY=cat(2,citY,citoY(k,:));
citC=cat(2,citC,citoC(k,:));
citR=cat(2,citR,citoR(k,:));
end


clear tgl frt fy fc tgl2 tgl3
n=64;
frt=Mfret(n,:);
fy=ffluoy(n,:);
fc=ffluoc(n,:);
fr=ffluor(n,:);
for k=n:n+4
frt=cat(2,frt,Mfret(k,:));
fy=cat(2,fy,ffluoy(k,:));
fc=cat(2,fc,ffluoc(k,:));
fr=cat(2,fr,ffluor(k,:));
end
for i=size(fr,2):-1:1
    if fr(i)==0
        fr(i)=[];
        fc(i)=[];
        fy(i)=[];
    end
end
% figure, title(['CFP vs FRET in frames ' num2str(n) ' to ' num2str(n+4) ' File 2_ 001']), hold on, scatter(fc(:),frt(:), 'b','.'), axis([60 170 0 1]), xlabel('CFP'), ylabel('FRET')
% figure, title(['YFP vs FRET in frames ' num2str(n) ' to ' num2str(n+4) ' File 2_ 001']), hold on, scatter(fy(:),frt(:), 'b','.'), axis([60 170 0 1]), xlabel('YFP'), ylabel('FRET')
figure, title(['Attivita` FRET (YFP) vs YFP (magenta) in frames ' num2str(n) ' to ' num2str(n+4) ' File 2_ 001']), hold on, scatter(fr(:),fy(:), 'b','.'), axis([0 60 0 110]), xlabel('magenta'), ylabel('YFP'),
hold on, scatter(citR(:),citY(:), 'r','.'), hold off
figure, title(['Attivita` FRET (YFP) vs CFP (magenta) in frames ' num2str(n) ' to ' num2str(n+4) ' File 2_ 001']), hold on, scatter(fc(:),fy(:), 'b','.'), axis([20 220 0 110]), xlabel('CFP'), ylabel('YFP'),
hold on, scatter(citC(:),citY(:), 'r','.'), hold off


% h = get(0,'children');
% for i=1:length(h)
%   saveas(h(i), ['sonda+vps34in1+trf_2_002_' num2str(i)], 'svg');
% end