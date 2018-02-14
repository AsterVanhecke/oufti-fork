function sgn = getOneSignalMedian(mesh,mbox,img)
    % This function returns the median signal within each segment of the
    % mesh. Based on getOneSignalMedian.
    % 
    
    if size(mesh,1)>1
        img2 = imcrop(img,mbox);
        sgn=nan(size(mesh,1)-1 , 1);
        for ii=1:size(mesh,1)-1
            % Create a polygon representing one segment
            plgx = [mesh(ii,[1 3]) mesh(ii+1,[3 1])]-mbox(1)+1;
            plgy = [mesh(ii,[2 4]) mesh(ii+1,[4 2])]-mbox(2)+1;
            % Make a pixel mask from the segment polygon.
            mask = poly2mask(plgx,plgy,size(img2,1),size(img2,2));
            % Count the number of pixels in the mask.
            s = sum(sum(mask));
            % Signal: Median of pixel intensities within the segment.
            if s>0, s=nanmedian(img2(mask)); end % *polyarea(plgx,plgy)/s
            sgn(ii) = s;
        end
    else
        sgn = [];
    end
end