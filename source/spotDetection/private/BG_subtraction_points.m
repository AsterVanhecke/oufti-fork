function out=BG_subtraction_points(cellData,rawImage,bgr,sigma,varargin)

mesh=cellData.mesh;
mesh(:,1)=mesh(:,1)-cellData.box(1)+1;% + align(1);
mesh(:,2)=mesh(:,2)-cellData.box(2)+1;% + align(2);
mesh(:,3)=mesh(:,3)-cellData.box(1)+1;% + align(1);
mesh(:,4)=mesh(:,4)-cellData.box(2)+1;% + align(2);

rawImage_msr_Intensity=rawImage;
if ~isempty(varargin)
    maskOfNewImage=varargin{1};
    rawImage_msr_Intensity=double(rawImage);
    rawImage_msr_Intensity(maskOfNewImage)=NaN;
end

% bgr=median(median(image));
% useMedian=false;
% if useMedian
%     signalArray = getOneSignalMedian(double(cellData.mesh),double(cellData.box),double(image)); 
%     segIntensity=signalArray-double(bgr);
%     segIntensity(segIntensity<0)=0;
%     segIntensity=segIntensity./1.73*2*pi;
% else
signalArrayS = getOneSignalM(double(mesh),[0 0 (cellData.box(3:4)+1)],double(rawImage_msr_Intensity)-double(bgr),1);
% end

diffX = mesh(:,4) - mesh(:,2);
diffY = mesh(:,3) - mesh(:,1);
widths = sqrt(diffX .^ 2 + diffY .^ 2);
segWidths=widths(1:end-1)+diff(widths)./2;
% if ~useMedian
segVols=pi*(segWidths/2).^2;
IperV= nanmedian(signalArrayS./segVols); % Calculate median of Intensity per Volume of each segment
% IperV= sum(signalArrayS)./sum(segVols);% calculate mean Intensity per Volume
if sum(~isnan(signalArrayS))<1 % If signalArrayS consists solely of NaN's. (or less than 3?)
    signalArrayS = getOneSignalM(double(mesh),[0 0 (cellData.box(3:4)+1)],double(rawImage-bgr),1);
    IperV= nanmedian(signalArrayS./segVols);
    warning('Too many spots to measure background Intensity')
end
segIntensity=IperV*segVols; % Total segment intensity is proportional to Volume
% end

cLine=([mean([mesh(:,1), mesh(:,3)],2) , mean([mesh(:,2), mesh(:,4)],2) ]);
segCent=cLine(1:end-1,:)+diff(cLine)./2;

% Vector defining width axis in a segment:
% center of a segment edge (between two consecutive points on the mesh) minus segment center.
segWvect=mesh(1:end-1,1:2)+diff(mesh(:,1:2))./2 - segCent;

segLengthVec =diff(cLine);

meshDiff=diff(mesh);
l1=sqrt(meshDiff(:,1).^2 + meshDiff(:,2).^2);
l2=sqrt(meshDiff(:,3).^2 + meshDiff(:,4).^2);
% normalize l1 and l2
lboth=[l1,l2];
l1= lboth(:,1)./sum(lboth,2)*2; l2= lboth(:,2)./sum(lboth,2)*2;

% Constants sampling density:
samplingD=5;
pointsPerI=10;

Nlines = ceil(segWidths .* samplingD); % ceil ? Always create at least one point
% N is 42 by 1
maxPtsGuess=round(length(segWidths)* max(segIntensity) *pointsPerI);
allX=nan(1,maxPtsGuess);
allY=nan(1,maxPtsGuess);
counter=1;
for segIdx=1:length(segWidths)
    % d can be calculated in a cell array:
    % is now 1 by N
    d=linspace(-1,1, Nlines(segIdx) ); % create X sampling for cross-section: relative position along the diameter of the cell
    y=real(2*sqrt(1^2-d.^2)); % d could be called from cell array % 1 by N
    y=y/sum(y); % Normalize y. This way total intensity is distributed throughout segment.
    % X and Y coordinates: of  multiply each d with the segment Width vector and add to the segment center.
    Xd=d*segWvect(segIdx,1)+segCent(segIdx,1); % 1 by N
    Yd=d*segWvect(segIdx,2)+segCent(segIdx,2); % 1 by N
    % correct for segments being longer at one edge due to cell curvature.
    segLenCorr=linspace( l2(segIdx) , l1(segIdx) ,length(d));
    % at each point, create regularly spaced point in the direction of the length axis.
    % Correct for curvature: segments don't have uniform size along the length dimension. (they are not rectangular)
    % Sampling: start at -0.5+ 1 point, in order not to have overlapping (double) sampling at segment edges 
    Nsampling=round(y.*segLenCorr.* segIntensity(segIdx).*pointsPerI);
    for dIdx=1:length(d)
        if Nsampling(dIdx)>0
            newX= Xd(dIdx) + segLenCorr(dIdx).* linspace(-0.5+1/Nsampling(dIdx),0.5, Nsampling(dIdx)) *segLengthVec(segIdx,1);
            newY= Yd(dIdx) + segLenCorr(dIdx).* linspace(-0.5+1/Nsampling(dIdx),0.5, Nsampling(dIdx)) *segLengthVec(segIdx,2);
            allX(counter:counter+Nsampling(dIdx)-1)= newX;
            allY(counter:counter+Nsampling(dIdx)-1)= newY;
            counter=counter+Nsampling(dIdx);
        end
    end
end

Xedges = (0:size(rawImage,2)) +0.5; 
Yedges = (0:size(rawImage,1)) +0.5;
Values=histcounts2(allX,allY,Xedges,Yedges);

% % median y ~1.73:
% d=linspace(-1,1, 10^6);
% y=real(2*sqrt(1^2-d.^2));
% median(y)

% The units of space are pixels here.
RescaleFactor=1/(pointsPerI);

% PSFsigma=PSFwidth/2;
% pxlSize=0.0636;
% sigma=PSFsigma./pxlSize;
Background = imgaussfilt(Values',sigma); % Blur, preserves total intensity

out.backgroundImage=RescaleFactor*Background;
out.height=Values';
out.bg_subtracted=double(rawImage)-out.backgroundImage;
out.RescaleFactor=RescaleFactor;

% %plot inline -s 1800,600
% figure, 
% subplot(1,3,1)
% imagesc(rawImage), axis image
% colorbar, title('Raw image')
% subplot(1,3,2)
% imagesc(RescaleFactor*Background), axis image
% colorbar, title('Background')
% subplot(1,3,3)
% imagesc(double(rawImage)-RescaleFactor*Background), axis image
% colorbar, title('Raw image after background subtraction')