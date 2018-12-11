function  [intSlo, stat, Pows, Deviants,  stat0, intSlo0] = fitPowerLaw3steps(XX,YY, robRegMeth,  doPlot, thisCol)
%% 
% fit a line on log-freq vs log-PSD, excluding peaks and their base.
% USAGE EXAMPLE: 
% % first compute the PSD
% epLen= 2* sRate; epShift= 1*srate;numFFT=[];
%  [myPSD,frex]= pwelch( myEEGch  , epLen, epShift,numFFT, sRate); 
%  frBand=[1 40];
%  frBins= dsearchn( frex, frBand(1) ):  dsearchn( frex, frBand(2));
% XX= frex(frBins);
% YY= myPSD(frBins);
% robRegMeth= 'ols';
% doPlot= 1;
% thisCol= [0 0 1];
%  [intSlo, stat, Pows, Deviants,  stat0, intSlo0] = fitPowerLaw3steps(XX,YY, robRegMeth,  doPlot, thisCol)
%%%%%%%%% INPUT
% XX ---> are the frequency bins of interest (hz)
% YY ---> are the PSD relative to the bins of interest (mV^2 /Hz)
% robRegMeth --->is the method for robust regression fit (default: ols)
%                type help robustfit for more details
% doPlot  --->    0 (default) means no plot ; 
%                 1 is plot XX and YY on log-log scale;  
%                 2 is plot log(XX) and log(YY) on linear scale; 
% thisCol --->    if plotting, is the rgb triplet
%%%%%%%%% OUTPUT
% intSlo(1) intercept of 2nd (final) powerLaw Fit
% intSlo(2) slope of 2nd (final) powerLaw Fit  i.e. the spectral exponent

% stat is a structure with informations on the 2nd fit, (see 2nd argument of robustfit)

% Pows is a structure with
%           pred: 1*length(frBins) doubles, predicted PSD by the power-law fit
%            obs: 1*length(frBins) doubles, observed PSD    = YY
%            res: 1*length(frBins) doubles, residuals (pred-obs) PSD
%           frex: 1*length(frBins) doubles, observed frex   = XX
%       meanPred: 1 double, mean of the predicted PSD
%        meanObs: 1 double, mean of the observed PSD
%        meanRes: 1 double, mean of the residuals of the PSD
%     meanResPos: 1 double, mean of the residuals of the PSD (exclude negative values)
%       cenFrObs: 1 double, central Frequency of the observed PSD
%       cenFrRes: 1 double, central Frequency of the residual PSD

% Deviants is a structure with 
%       res2: 1*length(frBins) doubles, 
%        res: 1*length(frBins) doubles, 
%        rej:  length(frBins)*1 logicals, 1 if bins are rejected, 0 if kept   for fitting the 2nd power-law
%         fr: min and max freqs considered
%       thre: threshold used for bins adjacent to the peaks
%     clusts: structure with information on the contiguous freq. bins 

% stat0 is a structure with informations on the 1st fit, (see 2nd argument of robustfit)

% intSlo(1) intercept of 2nd (final) powerLaw Fit
% intSlo(2) slope of 2nd (final) powerLaw Fit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HOW DOES IT WORK?
% transform in log-log and upsample by 4 times
% fit a 1st line, find all residual >0
% small peaks (whose residuals are smaller than 1*mad(residuals)) do not count as peaks
% find clusters of rejected frequency bins 
% consider only clusters of bins where there is a peak
% i.e. reject large enough peaks and their base (adjacent residuals >0)
% fit a second line, on the remaining residuals (those closer to a power-law)

if exist('robRegMeth','var')
else
    robRegMeth= 'ols';
end
%% LOG TRANSFORM
%%%%%%%%%%%%%%%%%
% vectorize frex and psd to avoid dimension mismatch later
XX= XX(:)';
YY= YY(:)';

 X= log10(XX);
 Y= log10(YY);

 %%% INTERPOLATE IN LOG-LOG: X, Y --> Xi, Yi
 XXi= logspace( X(1), X(end),  length(X)*4); % 2^10*2;--> auto-chosen
  Xi= log10(XXi); 
Yi= interpn( X, Y, Xi   );% ,'spline'
 YYi= 10.^(Yi); 

%% STEP 1, FIT 1st LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MODEL 1: minimize Y residuals. OLS +variants
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ intSlo0, stat0] = robustfit(Xi,Yi,robRegMeth);
    YRes0= stat0.resid;
    YPred0=  ( intSlo0(1)+intSlo0(2)*(Xi))';

%% FIND DEVIANT RESIDUALS
%%%%%%%%%%%%%%%%%%%%%%%%%
threRes=0; 
boolYdev= YRes0 > threRes ;% 
%% ONLY CONSIDER DEVIANTS WHERE ALSO POWER-PEAKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      %(I.E. if power is high, but it decays so that it does not peak (e.g.peak is lower than earliest freq bin) --> not deviant
% peaks are searched on log-power .... (why not on residuals? avoid small bumps on res)
[pks0,locs0] = findpeaks( Yi ); %, 'MINPEAKDISTANCE', 8, 'NPEAKS', 10, 'SORTSTR','descend') ;

% REJECT SMALL PEAKS
% peaks are rejected if their residuals (detrended log-power) are low ... (why not on log-power? there is no meaningful threshold on log-pow)
pks=pks0; locs=locs0;

threResPks= 1* mad(YRes0,1); %    median(YRes0) + 1 *mad(myres,1); 
threResPks= max(threResPks, .1);

rejectPks= (YRes0(locs) < threResPks )';
pks(rejectPks)=[];
locs(rejectPks)=[];

% now exclude deviant-residuals that do not contain any peak
Clusts = bwconncomp(boolYdev);
includeClust= false( 1, Clusts.NumObjects);
for cl= 1: Clusts.NumObjects
    thisIdx=Clusts.PixelIdxList{cl};
    thisclust= Xi( thisIdx );
    thisloc= (Xi(locs));
    includeClust(cl)=   sum( ismember ( thisclust , thisloc));
    boolYdev( thisIdx )=includeClust(cl);
end

%% STEP 2 ,2nd line after excluding peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % MODEL 1: minimize Y residuals. OLS +variants
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try % try an OLS regression first
        [ intSlo, stat] = robustfit(Xi(~boolYdev),Yi(~boolYdev),robRegMeth);
    catch % if not enough points for regression,(too small freq.band or deviations with too much band-width) just do not exclude any peak ans do ols
        [ intSlo, stat] = robustfit(Xi,Yi,'ols');
    end
    %%%% GET RESIDUALS ALSO AT EXCLUDED POINTS
    YPred0= (Xi.*intSlo(2) +intSlo(1));
    YRes2= Yi - YPred0;    %     YRes2= stat.resid;


%% STORE
stat.frBins= [ XXi(1)   XXi(end)  ];
stat.myName= [ num2str(round( stat.frBins(1) )) '-'  num2str(round(  stat.frBins(2)) ) ' Hz: '   num2str(intSlo(2)) ]; %num2str(intSlo')
stat.robRegMeth=robRegMeth;

Clusts.includeClust=includeClust;
Clusts.pks=pks;
Clusts.locs=locs;
Clusts.threResPks=threResPks;
Clusts.pks0=pks0;
Clusts.locs0=locs0;
Clusts.rejectPks0=rejectPks;

Deviants.res2= YRes2;
Deviants.res= YRes0;
Deviants.rej= boolYdev;
Deviants.fr= stat.frBins;
Deviants.thre= threRes;
Deviants.clusts=Clusts;

YYpred= 10.^( intSlo(1)+intSlo(2)*(Xi))';
resPow= YYi- YYpred' ; % figure; plot(XXi, resPow)

Pows.pred= YYpred';
Pows.obs= YYi;
Pows.res= resPow;
Pows.frex= XXi;

Pows.meanPred= nanmean( YYpred);
Pows.meanObs= nanmean( YYi);
Pows.meanRes= mean( resPow  ); %YYpred

idxPos= ( resPow) >0;%YYpred
Pows.meanResPos= mean( 10.^(Yi(idxPos))- YYpred(idxPos)'); %YYpred(idxPos)'
Pows.cenFrObs= sum(XXi.*YYi)./sum(YYi);
Pows.cenFrRes= sum(XXi.*resPow)./sum(resPow);

%% PLOTS
%%%%%%%%%  
%   figure
if ~exist('doPlot','var')
    doPlot=0;
end
if doPlot  
%     figure;
    if ~exist('thisCol','var')
        thisCol= 'k';
    end
    myLW=2;
end
if doPlot ==1 % PLOT ON LOG-LOG
    
    loglog(XXi, YYi, ':','color', thisCol); hold on;
    
    xPlot= XXi; xPlot(boolYdev)= nan;
    loglog(xPlot, 10.^(Yi ), '-',   'LineWidth',myLW, 'color', thisCol); hold on;
    %%% PLOT FIRST SLOPE
    YYpred0= 10.^( intSlo0(1)+intSlo0(2)*(Xi))';

    %%% PLOT SECOND (FINAL) SLOPE
    hp(1)= plot( XXi([1 end]), YYpred([1 end]),'-','LineWidth',myLW,'color', thisCol ); 
     legend(hp,stat.myName);
    
    hp(2)= plot( XXi([1 end]), YYpred0([1 end]),'--','LineWidth',myLW,'color', [.5 .5 .5] ); 

      %%%% DOTS ON PEAKS 
%     xPlot= XXi; xPlot(~boolYdev)= nan;
%     loglog( xPlot, 10.^(Yi ), 'o','markersize',4, 'color', thisCol); hold on;
%     loglog(xPlot(Clusts.locs), 10.^(Yi(Clusts.locs) ) ,'.', 'color', thisCol); %.k

elseif doPlot==2 % PLOT ON LINEAR SCALES (but use log(X) log(Y) )
   
    plot( (Xi),  (Yi), ':','color', thisCol); hold on;
    
    xPlot=  (Xi); xPlot(boolYdev)= nan;
    plot(xPlot,  (Yi ), '-', 'linewidth',myLW,'color', thisCol);   %'-'
    xPlot=  (Xi); xPlot(~boolYdev)= nan;
    plot( xPlot,  (Yi ), ':', 'linewidth',2, 'color', thisCol ); % 'o' 'markersize',4,


    Ypred= ( intSlo(1)+intSlo(2)*(Xi));
     hp= plot(  (Xi), Ypred,'-','LineWidth',myLW,'color', thisCol );
     legend(hp,stat.myName);
     
     %%% PLOT FIRST SLOPE
     Ypred0=  ( intSlo0(1)+intSlo0(2)*(Xi))';
     hp(1)= plot( Xi([1 end]), Ypred0([1 end]),'--','LineWidth',myLW,'color', [.5 .5 .5] );
    
       %%%% DOTS ON PEAKS 
%  plot(xPlot(Clusts.locs), Yi(Clusts.locs) ,'.','color',thisCol);
%     axis tight;
    
end
%%%%%%%%%%%%%%
%% debug PLOTS
%%%%%%%%%%%%%%
debugPlot=0;
if debugPlot ==1
    figure;
    sp(1)=   subplot(121);  % loglog(XXi  , Yi ,':k') 
    loglog(XX  , YY ,':k');hold on;      loglog(XXi  , 10.^(YPred0) ,':b')
    plot (XXi(~boolYdev) ,  (YYi(~boolYdev) ),'.b');
    YYthrePks=   10.^( intSlo0(1)+threResPks+intSlo0(2)*(Xi))';
    loglog(XXi  ,YYthrePks ,':k');
    loglog(XXi  , (YYpred) ,':r')
    
    sp(2)=  subplot(122);     plot(XXi, YRes0,':k'); hold on;
    plot (XXi(~boolYdev) ,  (YRes0(~boolYdev) ),'.b'); axis tight;
    plot(xlim,[threRes threRes],':b');  plot(xlim,[threResPks threResPks],':k');
    plot(XXi  , YRes2 ,':r')
    
    try
        subplot(121);
        plot(XXi(locs0), YYi(locs0),'.','color',[.5 .5 .5],'markersize',20 );
        subplot(122);
        plot(XXi(locs0), YRes0(locs0),'.','color',[.5 .5 .5],'markersize',20);
        
        subplot(121);
        plot(XXi(locs), YYi(locs),'.k' ); plot(XXi(locs(1)),YYi(locs(1)),'.r','markersize',20),axis tight;
        subplot(122);
        plot(XXi(locs), YRes0(locs),'.k'); plot(XXi(locs(1)), YRes0(locs(1)),'.r','markersize',20)
        linkaxes(sp,'x');
    catch
    end
end
