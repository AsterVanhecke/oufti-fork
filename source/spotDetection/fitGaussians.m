function [spotStructure,dispStructure] = fitGaussians(cellData,numberOfSpots,rawImage,newRawImage,backgroundRawImage,params,dilatedCellContour)
% TODO: document

%% Initialize parameters
if isfield(params,'postMinHeight')
    postFitMinHeight = params.postMinHeight;
    postFitMinWidth  = params.postMinWidth;
    postFitMaxWidth  = params.postMaxWidth;
    postFitError     = params.postError;
else
    postFitMinHeight = 0.0;
    postFitMinWidth  = 0.5;
    postFitMaxWidth  = 10;
    postFitError     = 0.0;
end
maxRadius       = params.fitRadius;
multGauss       = params.multGauss;
I = NaN;
spotStructure.l                            = [];
spotStructure.d                            = [];
spotStructure.x                            = [];
spotStructure.y                            = [];
spotStructure.positions                    = [];
spotStructure.adj_Rsquared                 = [];
spotStructure.confidenceInterval_x_y       = [];
spotStructure.w                            = [];
spotStructure.h                            = [];
spotStructure.magnitude                    = [];
spotStructure.b                            = [];
dispStructure.adj_Rrsquared                = [];
dispStructure.w                            = [];
dispStructure.h                            = [];

[rows,columns] = size(rawImage);
[newRows, newColumns] = meshgrid(1:columns,1:rows);

%% Identify spots to fit
counter = 0;
tempRawImage = bwlabel(newRawImage);
pixelNumbers = cellfun(@numel,numberOfSpots.PixelIdxList); % Counts nr of pixels in each spot
indexToGoodSpots = pixelNumbers > 2; % remove spots with less than two pixels
meanIndex = mean(pixelNumbers(indexToGoodSpots));
pixelNumMeanStd = meanIndex+std(pixelNumbers(indexToGoodSpots))/2;
%backgroundRawImage = mean(mean(rawImage)); % Use from input instead

%% Start loop through spots
for k = 1:numberOfSpots.NumObjects
    if numel(numberOfSpots.PixelIdxList{k}) < params.minRegionSize
        continue;
    end
    try
        tempNewRawImage = k == tempRawImage; % Get mask with pixels in spot k.
        tempNewRawImage = tempNewRawImage.*im2double(rawImage); % Acquire raw pixel values
        indexToSpots = numberOfSpots.PixelIdxList{k}; % Acquire pixel indices
        % Get first guess for location of maximum
        [~,id] = max(rawImage(indexToSpots)); % find pixel in spot with highest value.
        peakValueOfSpots = indexToSpots(id);
        % peakValueOfSpots is the pixel index (linear indexing), the
        % following lines convert it to row/column.
        rowPosition = (rem(peakValueOfSpots-1,rows)+1);
        columnPosition = ceil(peakValueOfSpots./rows);
        % Estimate position of peak with subpixel precision
        positionXY = positionEstimate(rawImage(rowPosition-1:rowPosition+1,columnPosition-1:columnPosition+1)) + [columnPosition-1,rowPosition-1] -1;
        % Convert to row/column
        rowPosition = positionXY(2);
        columnPosition = positionXY(1);
        if  pixelNumbers(k) > pixelNumMeanStd || numberOfSpots.NumObjects == 1 % if spot has more pixels than most spots (mean+std) or is only spot
            % Convolute with ones(5), to expand mask
            tempMask = conv2(double(tempNewRawImage),ones(5),'same');
            tempMask = logical(tempMask);
        else
            tempMask = ((newRows-columnPosition).^2 + (newColumns-rowPosition).^2).^(1/2) <= maxRadius;            
        end
        indexToMask = find(tempMask ==1);
        tempNewRawImage =rawImage(tempMask); % mask raw Image with expanded mask.
        
        tempNewRawImage1 = uint16(tempMask.*double(rawImage));
        [yy,xx] = find(tempNewRawImage1); % Get row and columns for pixels in tempNewRawImage1
        points = [xx yy];
        
        %% Prepare parameters for fitoptions
        % for each pixel in the spot, calculate the distance to the other
        % pixels and return the largest distance, to estimate max Width.
        [d,~]  = pdist2(points,points,'euclidean','largest',1);
        maxWidthEstimate = ceil(max(d));
        indexPeakDistanceXY = indexToMask;
        % redo position estimate, on window around new rowPosition and
        % columnPosition (from previous position estimate)
        positionXY = positionEstimate(rawImage(rowPosition-1:rowPosition+1,columnPosition-1:columnPosition+1)) + [columnPosition-1,rowPosition-1] -1;
        % make 3 copies of position for the fit.
        positionX = positionXY(1);
        positionY = positionXY(2);
        positionX1 = positionX;
        positionY1 = positionY;
        positionX2 = positionX;
        positionY2 = positionY;
        positionX3 = positionX;
        positionY3 = positionY;
        
        if ~inpolygon(positionX,positionY,...
                dilatedCellContour(:,1),dilatedCellContour(:,2))
            continue;
        end
        
        % Estimate height as max pixel value within spot - background
        heightEstimate = (max(rawImage(indexPeakDistanceXY))-mean(mean(rawImage)));
        
        widthEstimate  = 1.5; % widthEstimate is hardcoded as 1.5?
        rowPosition = (rem(indexPeakDistanceXY-1,rows)+1);
        columnPosition = ceil(indexPeakDistanceXY./rows);
        
        % Make 3 copies of height and width. % Removed background
        heightEstimate1 = heightEstimate;
        widthEstimate1 = widthEstimate;
        heightEstimate2 = heightEstimate;
        widthEstimate2 = widthEstimate;
        heightEstimate3 = heightEstimate;
        widthEstimate3  = widthEstimate;
    catch
        continue;
    end
    
    %% Create fitoptions
    % Initialize fit options for multiple 2D-Gauss fits.
    % gauss2dFitOptionsN corresponds to fit with N Gaussians.
    gauss2dFitOptions1 = fitoptions('Method','NonlinearLeastSquares','Algorithm','Levenberg-Marquardt',...
        'Lower',[0,0,0],...
        'Upper',[Inf,Inf,maxWidthEstimate],'MaxIter', 600,...
        'Startpoint',[backgroundRawImage,heightEstimate,...
        widthEstimate1,positionX,positionY]);
    gauss2dFitOptions2 = fitoptions('Method','NonlinearLeastSquares','Algorithm','Levenberg-Marquardt',...
        'Lower',[0,0,0,-Inf,-Inf,0,0,0,-Inf,-Inf],...
        'Upper',[Inf,Inf,maxWidthEstimate,Inf,Inf,Inf,Inf,maxWidthEstimate,Inf,Inf],'MaxIter',600,...
        'Startpoint',[backgroundRawImage,heightEstimate,...
        widthEstimate,positionX,positionY,...
        heightEstimate1,...
        widthEstimate1,positionX1,positionY1]);
    
    gauss2dFitOptions3 = fitoptions('Method','NonlinearLeastSquares','Algorithm','Levenberg-Marquardt',...
        'Lower',[0,0,0,-Inf,-Inf,0,0,0,-Inf,-Inf,0,0,0,-Inf,-Inf],...
        'Upper',[Inf,Inf,maxWidthEstimate,Inf,Inf,Inf,Inf,maxWidthEstimate,Inf,Inf,...
        Inf,Inf,maxWidthEstimate,Inf,Inf],'MaxIter',600,...
        'Startpoint',[backgroundRawImage,heightEstimate,...
        widthEstimate,positionX,positionY,...
        heightEstimate1,...
        widthEstimate1,positionX1,positionY1,...
        heightEstimate2,...
        widthEstimate2,positionX2,positionY2]);
    
    gauss2dFitOptions4 = fitoptions('Method','NonlinearLeastSquares','Algorithm','Levenberg-Marquardt',...
        'Lower',[0,0,0,-Inf,-Inf,0,0,0,-Inf,-Inf,0,0,0,-Inf,-Inf,0,0,0,-Inf,-Inf],...
        'Upper',[Inf,Inf,maxWidthEstimate,Inf,Inf,Inf,Inf,maxWidthEstimate,Inf,Inf,...
        Inf,Inf,maxWidthEstimate,Inf,Inf,Inf,Inf,maxWidthEstimate,Inf,Inf],'MaxIter',600,...
        'Startpoint',[backgroundRawImage,heightEstimate,...
        widthEstimate,positionX,positionY,...
        heightEstimate1,...
        widthEstimate1,positionX1,positionY1,...
        heightEstimate2,...
        widthEstimate2,positionX2,positionY2,...
        heightEstimate3,...
        widthEstimate3,positionX3,positionY3]);
    
    %% Create fittype's
    gauss2d1 = fittype(@(backgroundRawImage,heightEstimate,widthEstimate,positionX,positionY,x,y) ...
        backgroundRawImage+heightEstimate*exp(-(x-positionX).^2 ...
        /(2*widthEstimate^2)-(y-positionY).^2/(2*widthEstimate^2)),...
        'independent', {'x', 'y'},'dependent', 'z','options',gauss2dFitOptions1);
    
    gauss2d2 = fittype(@(backgroundRawImage,heightEstimate,widthEstimate,positionX,positionY,...
        heightEstimate1,widthEstimate1,positionX1,positionY1,x,y) ...
        backgroundRawImage+heightEstimate*exp(-(x-positionX).^2 ...
        /(2*widthEstimate^2)-(y-positionY).^2/(2*widthEstimate^2)) + ...
        heightEstimate1*exp(-(x-positionX1).^2 ...
        /(2*widthEstimate1^2)-(y-positionY1).^2/(2*widthEstimate1^2)),...
        'independent', {'x', 'y'},'dependent', 'z','options',gauss2dFitOptions2);
    
    gauss2d3 = fittype(@(backgroundRawImage,heightEstimate,widthEstimate,positionX,positionY,...
        heightEstimate1,widthEstimate1,positionX1,positionY1,...
        heightEstimate2,widthEstimate2,positionX2,positionY2,x,y) ...
        backgroundRawImage+heightEstimate*exp(-(x-positionX).^2 ...
        /(2*widthEstimate^2)-(y-positionY).^2/(2*widthEstimate^2)) + ...
        heightEstimate1*exp(-(x-positionX1).^2 ...
        /(2*widthEstimate1^2)-(y-positionY1).^2/(2*widthEstimate1^2)) + ...
        heightEstimate2*exp(-(x-positionX2).^2 ...
        /(2*widthEstimate2^2)-(y-positionY2).^2/(2*widthEstimate2^2)),...
        'independent', {'x', 'y'},'dependent', 'z','options',gauss2dFitOptions3);
    
    gauss2d4 = fittype(@(backgroundRawImage,heightEstimate,widthEstimate,positionX,positionY,...
        heightEstimate1,widthEstimate1,positionX1,positionY1,...
        heightEstimate2,widthEstimate2,positionX2,positionY2,...
        heightEstimate3,widthEstimate3,positionX3,positionY3,x,y) ...
        backgroundRawImage+heightEstimate*exp(-(x-positionX).^2 ...
        /(2*widthEstimate^2)-(y-positionY).^2/(2*widthEstimate^2)) + ...
        heightEstimate1*exp(-(x-positionX1).^2 ...
        /(2*widthEstimate1^2)-(y-positionY1).^2/(2*widthEstimate1^2)) + ...
        heightEstimate2*exp(-(x-positionX2).^2 ...
        /(2*widthEstimate2^2)-(y-positionY2).^2/(2*widthEstimate2^2)) + ...
        heightEstimate3*exp(-(x-positionX3).^2 ...
        /(2*widthEstimate3^2)-(y-positionY3).^2/(2*widthEstimate3^2)),...
        'independent', {'x', 'y'},'dependent', 'z','options',gauss2dFitOptions4);
    try
        MultipleGauss = {gauss2d1,gauss2d2,gauss2d3,gauss2d4};
    catch
    end
    
    %% Perform initial fit(s)
    try
        sfit = cell(1,4);
        gof  = cell(1,4);
        try
            if multGauss ~= 0
                for ii = 1:4
                    [sfit{ii},gof{ii}] = fit([columnPosition,rowPosition],double(tempNewRawImage),MultipleGauss{ii});
                    
                end
            else
                for ii = 1:1
                    [sfit{ii},gof{ii}] = fit([columnPosition,rowPosition],double(tempNewRawImage),MultipleGauss{ii});
                    
                end
            end
        catch
        end
        
        %% If there are many pixels, fit multiple Gaussians, starting from local maxima
        % If there are more pixels in the spot than can be expected from a single Gaussian,
        % and fitting multiple Gaussians is enabled.
        if (pixelNumbers(k) >= maxRadius^2*pi) && (multGauss == 1)
            try
                maxValues = findLocalMaximas(tempNewRawImage1);
                maxValues = unique(maxValues,'rows');
                num = size(maxValues,1); % Count nr of unique maxima
                % Perform fits with 1 Gaussian per found local maximum,
                % when the number of found maxima is 2, 3 or 4.
                switch num
                    case 1
                        % Do nothing
                    case 2
                        positionX = maxValues(1,1); % Use positions of local maxima
                        positionY = maxValues(1,2);
                        positionX1 = maxValues(2,1);
                        positionY1 = maxValues(2,2);
                        gauss2dFitOptions2.StartPoint(4) = positionX;
                        gauss2dFitOptions2.StartPoint(5) = positionY;
                        gauss2dFitOptions2.StartPoint(8) = positionX1;
                        gauss2dFitOptions2.StartPoint(9) = positionY1;
                        gauss2d2 = fittype(@(backgroundRawImage,heightEstimate,widthEstimate,positionX,positionY,...
                            heightEstimate1,widthEstimate1,positionX1,positionY1,x,y) ...
                            backgroundRawImage+heightEstimate*exp(-(x-positionX).^2 ...
                            /(2*widthEstimate^2)-(y-positionY).^2/(2*widthEstimate^2)) + ...
                            heightEstimate1*exp(-(x-positionX1).^2 ...
                            /(2*widthEstimate1^2)-(y-positionY1).^2/(2*widthEstimate1^2)),...
                            'independent', {'x', 'y'},'dependent', 'z','options',gauss2dFitOptions2);
                    case 3
                        positionX = maxValues(1,1);
                        positionY = maxValues(1,2);
                        positionX1 = maxValues(2,1);
                        positionY1 = maxValues(2,2);
                        positionX2 = maxValues(3,1);
                        positionY2 = maxValues(3,2);
                        gauss2dFitOptions3.StartPoint(4) = positionX;
                        gauss2dFitOptions3.StartPoint(5) = positionY;
                        gauss2dFitOptions3.StartPoint(8) = positionX1;
                        gauss2dFitOptions3.StartPoint(9) = positionY1;
                        gauss2dFitOptions3.StartPoint(12) = positionX2;
                        gauss2dFitOptions3.StartPoint(13) = positionY2;
                        
                        gauss2d3 = fittype(@(backgroundRawImage,heightEstimate,widthEstimate,positionX,positionY,...
                            heightEstimate1,widthEstimate1,positionX1,positionY1,...
                            heightEstimate2,widthEstimate2,positionX2,positionY2,x,y) ...
                            backgroundRawImage+heightEstimate*exp(-(x-positionX).^2 ...
                            /(2*widthEstimate^2)-(y-positionY).^2/(2*widthEstimate^2)) + ...
                            heightEstimate1*exp(-(x-positionX1).^2 ...
                            /(2*widthEstimate1^2)-(y-positionY1).^2/(2*widthEstimate1^2)) + ...
                            heightEstimate2*exp(-(x-positionX2).^2 ...
                            /(2*widthEstimate2^2)-(y-positionY2).^2/(2*widthEstimate2^2)),...
                            'independent', {'x', 'y'},'dependent', 'z','options',gauss2dFitOptions3);
                        gauss2dFitOptions2.StartPoint(4) = maxValues(1,1);
                        gauss2dFitOptions2.StartPoint(5) = maxValues(1,2);
                        gauss2dFitOptions2.StartPoint(8) = maxValues(2,1);
                        gauss2dFitOptions2.StartPoint(9) = maxValues(2,2);
                        gauss2d2 = fittype(@(backgroundRawImage,heightEstimate,widthEstimate,positionX,positionY,...
                            heightEstimate1,widthEstimate1,positionX1,positionY1,x,y) ...
                            backgroundRawImage+heightEstimate*exp(-(x-positionX).^2 ...
                            /(2*widthEstimate^2)-(y-positionY).^2/(2*widthEstimate^2)) + ...
                            heightEstimate1*exp(-(x-positionX1).^2 ...
                            /(2*widthEstimate1^2)-(y-positionY1).^2/(2*widthEstimate1^2)),...
                            'independent', {'x', 'y'},'dependent', 'z','options',gauss2dFitOptions2);
                        
                    case 4
                        gauss2dFitOptions2.StartPoint(4) = maxValues(1,1);
                        gauss2dFitOptions2.StartPoint(5) = maxValues(1,2);
                        gauss2dFitOptions2.StartPoint(8) = maxValues(2,1);
                        gauss2dFitOptions2.StartPoint(9) = maxValues(2,2);
                        
                        gauss2d2 = fittype(@(backgroundRawImage,heightEstimate,widthEstimate,positionX,positionY,...
                            heightEstimate1,widthEstimate1,positionX1,positionY1,x,y) ...
                            backgroundRawImage+heightEstimate*exp(-(x-positionX).^2 ...
                            /(2*widthEstimate^2)-(y-positionY).^2/(2*widthEstimate^2)) + ...
                            heightEstimate1*exp(-(x-positionX1).^2 ...
                            /(2*widthEstimate1^2)-(y-positionY1).^2/(2*widthEstimate1^2)),...
                            'independent', {'x', 'y'},'dependent', 'z','options',gauss2dFitOptions2);
                        positionX = maxValues(1,1);
                        positionY = maxValues(1,2);
                        positionX1 = maxValues(2,1);
                        positionY1 = maxValues(2,2);
                        positionX2 = maxValues(3,1);
                        positionY2 = maxValues(3,2);
                        positionX3 = maxValues(4,1);
                        positionY3 = maxValues(4,2);
                        gauss2dFitOptions3.StartPoint(4) = positionX;
                        gauss2dFitOptions3.StartPoint(5) = positionY;
                        gauss2dFitOptions3.StartPoint(8) = positionX1;
                        gauss2dFitOptions3.StartPoint(9) = positionY1;
                        gauss2dFitOptions3.StartPoint(12) = positionX2;
                        gauss2dFitOptions3.StartPoint(13) = positionY2;
                        
                        gauss2d3 = fittype(@(backgroundRawImage,heightEstimate,widthEstimate,positionX,positionY,...
                            heightEstimate1,widthEstimate1,positionX1,positionY1,...
                            heightEstimate2,widthEstimate2,positionX2,positionY2,x,y) ...
                            backgroundRawImage+heightEstimate*exp(-(x-positionX).^2 ...
                            /(2*widthEstimate^2)-(y-positionY).^2/(2*widthEstimate^2)) + ...
                            heightEstimate1*exp(-(x-positionX1).^2 ...
                            /(2*widthEstimate1^2)-(y-positionY1).^2/(2*widthEstimate1^2)) + ...
                            heightEstimate2*exp(-(x-positionX2).^2 ...
                            /(2*widthEstimate2^2)-(y-positionY2).^2/(2*widthEstimate2^2)),...
                            'independent', {'x', 'y'},'dependent', 'z','options',gauss2dFitOptions3);
                        gauss2dFitOptions4.StartPoint(4) = positionX;
                        gauss2dFitOptions4.StartPoint(5) = positionY;
                        gauss2dFitOptions4.StartPoint(8) = positionX1;
                        gauss2dFitOptions4.StartPoint(9) = positionY1;
                        gauss2dFitOptions4.StartPoint(12) = positionX2;
                        gauss2dFitOptions4.StartPoint(13) = positionY2;
                        gauss2dFitOptions4.StartPoint(16) = positionX3;
                        gauss2dFitOptions4.StartPoint(17) = positionY3;
                        
                        gauss2d4 = fittype(@(backgroundRawImage,heightEstimate,widthEstimate,positionX,positionY,...
                            heightEstimate1,widthEstimate1,positionX1,positionY1,...
                            heightEstimate2,widthEstimate2,positionX2,positionY2,...
                            heightEstimate3,widthEstimate3,positionX3,positionY3,x,y) ...
                            backgroundRawImage+heightEstimate*exp(-(x-positionX).^2 ...
                            /(2*widthEstimate^2)-(y-positionY).^2/(2*widthEstimate^2)) + ...
                            heightEstimate1*exp(-(x-positionX1).^2 ...
                            /(2*widthEstimate1^2)-(y-positionY1).^2/(2*widthEstimate1^2)) + ...
                            heightEstimate2*exp(-(x-positionX2).^2 ...
                            /(2*widthEstimate2^2)-(y-positionY2).^2/(2*widthEstimate2^2)) + ...
                            heightEstimate3*exp(-(x-positionX3).^2 ...
                            /(2*widthEstimate3^2)-(y-positionY3).^2/(2*widthEstimate3^2)),...
                            'independent', {'x', 'y'},'dependent', 'z','options',gauss2dFitOptions4);
                end
                
                MultipleGauss = {gauss2d1,gauss2d2,gauss2d3,gauss2d4};
                sfit = cell(1,4);
                gof  = cell(1,4);
                for ii = 1:4 %redo fits, overwriting old ones.
                    [sfit{ii},gof{ii}] = fit([columnPosition,rowPosition],double(tempNewRawImage),MultipleGauss{ii});
                end
            catch
            end
        end
        
        %% Decide how many Gaussians fit best
        if (pixelNumbers(k) <= maxRadius^2*pi) || (multGauss == 0)
            maxIndex = 1;
        else
            gofTemp = cell2mat(gof(~cellfun(@isempty,gof)));
            gofTemp = cat(1,gofTemp.adjrsquare); %get adjusted Rsquare for 1-4 Gaussian fits
            [~,maxIndex] = max(gofTemp); % find how many Gaussian fits gives the best adjusted Rsquare.
        end
        
        %% Assign output from the chosen fit. (in case of single Gaussian, refit first)
        sfitTemp = sfit{maxIndex};
        switch maxIndex
            case 1
                if sfitTemp.backgroundRawImage < 0 || gof{1}.adjrsquare < 0.3
                    gauss2dFitOptions1.weights = [];
                    gauss2d1 = fittype(@(backgroundRawImage,heightEstimate,widthEstimate,positionX,positionY,x,y) ...
                        backgroundRawImage+heightEstimate*exp(-(x-positionX).^2 ...
                        /(2*widthEstimate^2)-(y-positionY).^2/(2*widthEstimate^2)),...
                        'independent', {'x', 'y'},'dependent', 'z','options',gauss2dFitOptions1);
                    [sfitTemp,gof{1}] = fit([columnPosition,rowPosition],double(tempNewRawImage1(indexPeakDistanceXY)),gauss2d1);
                    
                end
                if (sfitTemp.heightEstimate/65535) < postFitMinHeight || ((sfitTemp.heightEstimate/65535) > 2 ) || ...
                        abs(sfitTemp.widthEstimate)*sqrt(2) < postFitMinWidth || abs(sfitTemp.widthEstimate)*sqrt(2) > postFitMaxWidth || ...
                        gof{1}.adjrsquare < postFitError || ~inpolygon(sfitTemp.positionX,sfitTemp.positionY,...
                        dilatedCellContour(:,1),dilatedCellContour(:,2))
                    continue;
                end
                counter = counter + 1;
                try
                    confidenceInterval = confint(sfitTemp);
                    confidenceInterval = confidenceInterval(:,4:5);
                catch
                    confidenceInterval = [];
                end
                
                xModelValue = cellData.box(1)+sfitTemp.positionX-1;
                yModelValue = cellData.box(2)+sfitTemp.positionY-1;
                if ~isfield(cellData,'steplength')
                    cellData = getextradata(cellData);
                end
                [l,d] = projectToMesh(cellData.box(1)-1+sfitTemp.positionX,...
                    cellData.box(2)-1+sfitTemp.positionY,cellData.mesh,cellData.steplength);
                for kk = 1:size(cellData.mesh,1)-1
                    pixelPeakX = [cellData.mesh(kk,[1 3]) cellData.mesh(kk+1,[3 1])] - cellData.box(1)+1;
                    pixelPeakY = [cellData.mesh(kk,[2 4]) cellData.mesh(kk+1,[4 2])] - cellData.box(2)+1;
                    if inpolygon(sfitTemp.positionX,sfitTemp.positionY,pixelPeakX,pixelPeakY)
                        I = kk;
                        break
                    end
                end
                
                Q = 2*abs(pi*sfitTemp.heightEstimate*sfitTemp.widthEstimate^2);
                spotStructure.l(counter) = l;
                dispStructure.w(counter) = sfitTemp.widthEstimate*sqrt(2);
                dispStructure.h(counter) = sfitTemp.heightEstimate; %/65535;
                spotStructure.w(counter) = sfitTemp.widthEstimate*sqrt(2); % Added width and height to spotStructure
                spotStructure.h(counter) = sfitTemp.heightEstimate;%/65535; %  divided by 2^16, 16bit rescale?
                spotStructure.magnitude(counter) = Q;
                spotStructure.b(counter)=sfitTemp.backgroundRawImage;
                spotStructure.d(counter) = d;
                spotStructure.x(counter) = xModelValue;
                spotStructure.y(counter) = yModelValue;
                spotStructure.positions(counter) = I;
                spotStructure.adj_Rsquared(counter) = gof{1}.adjrsquare;
                spotStructure.confidenceInterval_x_y{counter} = confidenceInterval;
                
            case 2
                b = [sfitTemp.backgroundRawImage sfitTemp.backgroundRawImage];
                w = [abs(sfitTemp.widthEstimate) abs(sfitTemp.widthEstimate1)];
                Q = [2*abs(pi*sfitTemp.heightEstimate*sfitTemp.widthEstimate^2) 2*abs(pi*sfitTemp.heightEstimate1*sfitTemp.widthEstimate1^2)];
                h = [sfitTemp.heightEstimate sfitTemp.heightEstimate1];
                posX = [sfitTemp.positionX sfitTemp.positionX1];
                posY = [sfitTemp.positionY sfitTemp.positionY1];
                try
                    confidenceInterval = confint(sfitTemp);
                    confidenceInterval = confidenceInterval(:,4:5);
                catch
                    confidenceInterval = [];
                end
                %------------------------------------------------------------------
                %if sfit.heightEstimate/65535 is less than heightCutoff(2) then do
                %not count spot as true spot and continue
                for jj = 1:2
                    if (h(jj)/65535) < postFitMinHeight || ((h(jj)/65535) > 1 ) || ...
                            w(jj)*sqrt(2) < postFitMinWidth || w(jj)*sqrt(2) > postFitMaxWidth || ...
                            gof{2}.adjrsquare < postFitError || ~inpolygon(posX(jj),posY(jj),...
                            dilatedCellContour(:,1),dilatedCellContour(:,2))
                        continue;
                    end
                    
                    counter = counter + 1;
                    xModelValue = cellData.box(1)-1+posX(jj);
                    yModelValue = cellData.box(2)-1+posY(jj);
                    
                    if ~isfield(cellData,'steplength')
                        cellData = getextradata(cellData);
                    end
                    [l,d] = projectToMesh(cellData.box(1)-1+posX(jj),...
                        cellData.box(2)-1+posY(jj),cellData.mesh,cellData.steplength);
                    
                    for kk = 1:size(cellData.mesh,1)-1
                        pixelPeakX = [cellData.mesh(kk,[1 3]) cellData.mesh(kk+1,[3 1])] - cellData.box(1)+1;
                        pixelPeakY = [cellData.mesh(kk,[2 4]) cellData.mesh(kk+1,[4 2])] - cellData.box(2)+1;
                        if inpolygon(posX(jj),posY(jj),pixelPeakX,pixelPeakY)
                            I = kk;
                            break
                        end
                    end
                    spotStructure.l(counter) = l;
                    dispStructure.w(counter) = w(jj)*sqrt(2);
                    dispStructure.h(counter) = h(jj); %/65535;
                    spotStructure.w(counter) = w(jj)*sqrt(2); % Added width and height to spotStructure
                    spotStructure.h(counter) = h(jj);%/65535; %  divided by 2^16, 16bit rescale?
                    spotStructure.magnitude(counter) = Q(jj);
                    spotStructure.b(counter)=b(jj);
                    spotStructure.d(counter) = d;
                    spotStructure.x(counter) = xModelValue;
                    spotStructure.y(counter) = yModelValue;
                    spotStructure.positions(counter) = I;
                    spotStructure.adj_Rsquared(counter) = gof{2}.adjrsquare;
                    spotStructure.confidenceInterval_x_y{counter} = confidenceInterval;
                    
                    
                end
            case 3
                b = [sfitTemp.backgroundRawImage sfitTemp.backgroundRawImage sfitTemp.backgroundRawImage];
                w = [abs(sfitTemp.widthEstimate) abs(sfitTemp.widthEstimate1) abs(sfitTemp.widthEstimate2)];
                Q = [2*abs(pi*sfitTemp.heightEstimate*sfitTemp.widthEstimate^2)  ...
                    2*abs(pi*sfitTemp.heightEstimate1*sfitTemp.widthEstimate1^2)...
                    2*abs(pi*sfitTemp.heightEstimate2*sfitTemp.widthEstimate2^2)];
                h = [sfitTemp.heightEstimate sfitTemp.heightEstimate1 sfitTemp.heightEstimate2];
                posX = [sfitTemp.positionX sfitTemp.positionX1 sfitTemp.positionX2];
                posY = [sfitTemp.positionY sfitTemp.positionY1 sfitTemp.positionY2];
                try
                    confidenceInterval = confint(sfitTemp);
                    confidenceInterval = confidenceInterval(:,4:5);
                catch
                    confidenceInterval = [];
                end
                %------------------------------------------------------------------
                %if sfit.heightEstimate/65535 is less than heightCutoff(2) then do
                %not count spot as true spot and continue
                for jj = 1:3
                    if (h(jj)/65535) < postFitMinHeight || ((h(jj)/65535) > 1 ) || ...
                            w(jj)*sqrt(2) < postFitMinWidth || w(jj)*sqrt(2) > postFitMaxWidth || ...
                            gof{3}.adjrsquare < postFitError || ~inpolygon(posX(jj),posY(jj),...
                            dilatedCellContour(:,1),dilatedCellContour(:,2))
                        continue;
                    end
                    
                    counter = counter + 1;
                    xModelValue = cellData.box(1)-1+posX(jj);
                    yModelValue = cellData.box(2)-1+posY(jj);
                    
                    if ~isfield(cellData,'steplength')
                        cellData = getextradata(cellData);
                    end
                    [l,d] = projectToMesh(cellData.box(1)-1+posX(jj),...
                        cellData.box(2)-1+posY(jj),cellData.mesh,cellData.steplength); %#ok<AGROW>
                    for kk = 1:size(cellData.mesh,1)-1
                        pixelPeakX = [cellData.mesh(kk,[1 3]) cellData.mesh(kk+1,[3 1])] - cellData.box(1)+1;
                        pixelPeakY = [cellData.mesh(kk,[2 4]) cellData.mesh(kk+1,[4 2])] - cellData.box(2)+1;
                        if inpolygon(posX(jj),posY(jj),pixelPeakX,pixelPeakY)
                            I = kk;
                            break
                        end
                    end
                    spotStructure.l(counter) = l;
                    dispStructure.w(counter) = w(jj)*sqrt(2);
                    dispStructure.h(counter) = h(jj); %/65535;
                    spotStructure.w(counter) = w(jj)*sqrt(2); % Added width and height to spotStructure
                    spotStructure.h(counter) = h(jj);%/65535; %  divided by 2^16, 16bit rescale?
                    spotStructure.magnitude(counter) = Q(jj);
                    spotStructure.b(counter)=b(jj);
                    spotStructure.d(counter) = d;
                    spotStructure.x(counter) = xModelValue;
                    spotStructure.y(counter) = yModelValue;
                    spotStructure.positions(counter) = d;
                    spotStructure.adj_Rsquared(counter) = gof{3}.adjrsquare;
                    spotStructure.confidenceInterval_x_y{counter} = confidenceInterval;
                    
                end
            case 4
                b = [sfitTemp.backgroundRawImage sfitTemp.backgroundRawImage sfitTemp.backgroundRawImage sfitTemp.backgroundRawImage];
                w = [abs(sfitTemp.widthEstimate) abs(sfitTemp.widthEstimate1) abs(sfitTemp.widthEstimate2) abs(sfitTemp.widthEstimate2)];
                Q = [2*abs(pi*sfitTemp.heightEstimate*sfitTemp.widthEstimate^2)  ...
                    2*abs(pi*sfitTemp.heightEstimate1*sfitTemp.widthEstimate1^2)...
                    2*abs(pi*sfitTemp.heightEstimate2*sfitTemp.widthEstimate2^2)...
                    2*abs(pi*sfitTemp.heightEstimate3*sfitTemp.widthEstimate3^2)];
                h = [sfitTemp.heightEstimate sfitTemp.heightEstimate1 sfitTemp.heightEstimate2 sfitTemp.heightEstimate3];
                posX = [sfitTemp.positionX sfitTemp.positionX1 sfitTemp.positionX2 sfitTemp.positionX3];
                posY = [sfitTemp.positionY sfitTemp.positionY1 sfitTemp.positionY2 sfitTemp.positionY3];
                try
                    confidenceInterval = confint(sfitTemp);
                    confidenceInterval = confidenceInterval(:,4:5);
                catch
                    confidenceInterval = [];
                end
                %------------------------------------------------------------------
                %if sfit.heightEstimate/65535 is less than heightCutoff(2) then do
                %not count spot as true spot and continue
                for jj = 1:4
                    if (h(jj)/65535) < postFitMinHeight || ((h(jj)/65535) > 1 ) || ...
                            w(jj)*sqrt(2) < postFitMinWidth || w(jj)*sqrt(2) > postFitMaxWidth || ...
                            gof{4}.adjrsquare < postFitError || ~inpolygon(posX(jj),posY(jj),...
                            dilatedCellContour(:,1),dilatedCellContour(:,2))
                        continue;
                    end
                    
                    counter = counter + 1;
                    xModelValue = cellData.box(1)-1+posX(jj);
                    yModelValue = cellData.box(2)-1+posY(jj);
                    if ~isfield(cellData,'steplength')
                        cellData = getextradata(cellData);
                    end
                    [l,d] = projectToMesh(cellData.box(1)-1+posX(jj),...
                        cellData.box(2)-1+posY(jj),cellData.mesh,cellData.steplength);
                    for kk = 1:size(cellData.mesh,1)-1
                        pixelPeakX = [cellData.mesh(kk,[1 3]) cellData.mesh(kk+1,[3 1])] - cellData.box(1)+1;
                        pixelPeakY = [cellData.mesh(kk,[2 4]) cellData.mesh(kk+1,[4 2])] - cellData.box(2)+1;
                        if inpolygon(posX(jj),posY(jj),pixelPeakX,pixelPeakY)
                            I = kk;
                            break
                        end
                    end
                    spotStructure.l(counter) = l;
                    dispStructure.w(counter) = w(jj)*sqrt(2);
                    dispStructure.h(counter) = h(jj); %/65535; %  divided by 2^16, 16bit rescale?
                    spotStructure.w(counter) = w(jj)*sqrt(2); % Added width and height to spotStructure
                    spotStructure.h(counter) = h(jj);%/65535; %  divided by 2^16, 16bit rescale?
                    spotStructure.magnitude(counter) = Q(jj);
                    spotStructure.b(counter)=b(jj);
                    spotStructure.d(counter) = d;
                    spotStructure.x(counter) = xModelValue;
                    spotStructure.y(counter) = yModelValue;
                    spotStructure.positions(counter) = d;
                    spotStructure.adj_Rsquared(counter) = gof{4}.adjrsquare;
                    spotStructure.confidenceInterval_x_y{counter} = confidenceInterval;
                end
        end
        
    catch err
        disp(err);
    end
end % end for k=...
%% Sort output by l-coordinate
%--------------------------------------------------------------------------
%sort all fields of spotStructure according to the sorted "spotStructure.l"
%array.
[~,index] = sort(spotStructure.l);
spotStructure.l = reshape(spotStructure.l(index),1,[]);
spotStructure.d = reshape(spotStructure.d(index),1,[]);
spotStructure.x = reshape(spotStructure.x(index),1,[]);
spotStructure.y = reshape(spotStructure.y(index),1,[]);
spotStructure.positions = reshape(spotStructure.positions(index),1,[]);
spotStructure.adj_Rsquared = reshape(spotStructure.adj_Rsquared(index),1,[]);
dispStructure.adj_Rsquared  = reshape(spotStructure.adj_Rsquared(index),1,[]);
dispStructure.w        = reshape(dispStructure.w(index),1,[]);
dispStructure.h        = reshape(dispStructure.h(index),1,[]);
spotStructure.w        = reshape(spotStructure.w(index),1,[]);
spotStructure.h        = reshape(spotStructure.h(index),1,[]);
spotStructure.magnitude= reshape(spotStructure.magnitude(index),1,[]);
spotStructure.b        = reshape(spotStructure.b(index),1,[]);
%--------------------------------------------------------------------------