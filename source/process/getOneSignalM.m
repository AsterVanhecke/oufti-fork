function sgn = getOneSignalM(mesh,mbox,img,rsz)
    % This function integrates signal confined within each segment of the
    % mesh. Two versions are provided. This one is a faster but less
    % precise MATLAB-based approximation.
    % 
    
    if size(mesh,1)>1
        img2 = imcrop(img,mbox);
        if rsz>1
            img2 = imresize(img2,rsz);
        end
        sgn=nan(size(mesh,1)-1 , 1);
        for ii=1:size(mesh,1)-1
            % Create a polygon representing one segment
            plgx = rsz*([mesh(ii,[1 3]) mesh(ii+1,[3 1])]-mbox(1)+1);
            plgy = rsz*([mesh(ii,[2 4]) mesh(ii+1,[4 2])]-mbox(2)+1);
            % Make a pixel mask from the segment polygon.
            mask = poly2mask(plgx,plgy,size(img2,1),size(img2,2));
            % Count the number of pixels in the mask.
            s = sum(sum(mask));
            % Signal: sum of pixel intensities, corrected for pixel
            % area/segment polygon area. The correction counteracts
            % artifacts due to the pixelized nature of the intensity
            % signal. This way the total intensity of the pixels in the
            % cell is equal to the sum of the measured segment intensities
            % I don't know what the resize variable is for.
            if s>0, s=(rsz^(-2))*sum(sum(img2(mask)))*polyarea(plgx,plgy)/s; end
            sgn(ii) = s;
        end
    else
        sgn = [];
    end
end