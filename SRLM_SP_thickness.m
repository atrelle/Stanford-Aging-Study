function [srlm_out,sp_out] = thickness_sg_auto(blind_id,data_folder,anat_ext,mask_ext,show)
% Example Call:
% [srlm_out,sp_out] = thickness_sg_auto('0440_01212015_004','/share/awagner/AM/7T_ROIS/JEFF/Unblinded_Scans/440_01212015_004','.nii','.nii.gz','1')
% blind_id = 'blind033'
% data_folder '/share/awagner/AM/7T_ROIS/FSL_Analysis/current/work' note:
% no slash at end
% anat_ext, mask_ext = extension for image files: '.nii', '.nii.gz', '.img'
% show = 0 or 1. 0 = do not show work. 1 = show work.
% run sp_srlm_database_script before this
% -------------------------------------------------------------------------
% assumes size of x dim is 1024 (right is <512; left is > 512) - sag
% 23may2016

warning('off')

splineBreaks = 1;
splineOrder = 6;

% Getting input
anat_file = [data_folder '/' blind_id anat_ext];
srlm_file = [data_folder '/' blind_id '-mask' mask_ext]; % blind033-mask.img
sp_file   = [data_folder '/' blind_id 'ctx-mask' mask_ext]; % blind033ctx-mask.img

% SP,SRLM LOOP
for srlm = [0 1]; % (0 = sp, 1 = srlm)

    sampPoints = -10:10;
    xx = linspace(-10,10,10001);
    fig=0;
    clear fse;

    %clear randFse;
    fse = readFileNifti(anat_file);

    results = zeros(16,2); % changed by sag
    
    for left_right = 1:2; % added by sag

        clear srMask; % moved by sag
        clear sampFse; % moved by sag

        switch srlm
            case 0 % sp
                srMask = readFileNifti(sp_file);
            case 1 % srlm
                srMask = readFileNifti(srlm_file);
        end

        switch left_right % this switch added by sag
            case 1 % left
                srMask.data(1:512,:,:) = 0; % set right values to 0
            case 2 % right
                srMask.data(513:1024,:,:) = 0; % set left values to 0
        end

        cc = bwconncomp(srMask.data);

        for objNum=1:cc.NumObjects; % objNum=1:2 - edited by sag

            [x,y,z] = ind2sub(cc.ImageSize,cc.PixelIdxList{objNum});
            slices = unique(z);
            for sl=1:numel(slices)
                current_slice_number = slices(sl); % add by sag

                disp([srlm left_right current_slice_number]);

                curX = x(z==slices(sl));
                curY = y(z==slices(sl));
                curSliceInds = cc.PixelIdxList{objNum}(z==slices(sl));
                xFov = min(curX)-10:max(curX)+10;
                yFov = min(curY)-10:max(curY)+10;
                curSliceIm = fse.data(:,:,slices(sl));
                im = curSliceIm(xFov,yFov);

                [sortX,sortY] = cniSort2dPoints(curX, curY);
                t = cumsum([0;sqrt(diff(sortX(:)).^2 + diff(sortY(:)).^2)]); 
                sx = splinefit(t,sortX,splineBreaks,splineOrder); 
                sy = splinefit(t,sortY,splineBreaks,splineOrder);
                dt = 0.5;
                tt = t(1):dt:t(end); 
                splineX = ppval(sx, tt);
                splineY = ppval(sy, tt);

                if show > 0
                    fig = fig+1;
                    figure(fig); clf;

                    imagesc(xFov,yFov,im');
                    axis equal tight xy;
                    colormap(gray(256));
                    colorbar;
                    hold on;
                    plot(curX,curY,'.',splineX,splineY,'-');
                end
                sampFse = zeros(21,numel(curX));

                for(pt=2:numel(splineX))
                    p = polyfit([splineX(pt-1) splineX(pt)],[splineY(pt-1)   splineY(pt)],1);
                    m = -1/p(1);
                    if(abs(m)>1e6),      n = [0 1];
                    elseif(abs(m)<1e-6), n = [1 0];
                    else, n = [1 m]; end
                    % Make it unit length
                    n = n/norm(n);
                    sampX = mean([splineX(pt-1),splineX(pt)])+n(1)*sampPoints;
                    sampY = mean([splineY(pt-1),splineY(pt)])+n(2)*sampPoints;
                    sampFse(:,pt) = interp2(double(fse.data(:,:,slices(sl))'), sampX, sampY, 'linear');
                    sampFse(:,pt) = 0.5 * (sampFse(:,pt) + flipdim(sampFse(:,pt),1));
                    if show > 0
                        plot(sampX, sampY,'y-'); 
                    end
                end
                
                if show > 0
                    title(['Left/Right: ' num2str(left_right)  '. Slice: ' num2str(current_slice_number)]); % added by sag 12dec16
                end
                
                yyShift = spline(sampPoints,mean(sampFse,2),xx);
                yyShift1 = spline(sampPoints(2:21)-0.5,diff(mean(sampFse,2)),xx);
                if show > 0
                    fig=fig+1;
                    figure(fig); clf;
                    plot(sampPoints(:),mean(sampFse,2),'bo',xx,yyShift,'b-')
                    hold on
                    plot(sampPoints(2:21)-0.5,diff(mean(sampFse,2))+mean(yyShift),'ro',xx,yyShift1+mean(yyShift),'r-')
                end
                if srlm==1
                    [ngVal,ngInd] = min(yyShift1(3000:5000));
                    [poVal,poInd] = max(yyShift1(5001:7000));
                else
                    [ngVal,ngInd] = max(yyShift1(1000:5000));
                    [poVal,poInd] = min(yyShift1(5001:9000));
                end
%                 results(slices(sl),1+(median(curX)>(cc.ImageSize(1)/2)))=fse.pixdim(1)*20*(poInd-ngInd+4000-2000*srlm)/10001;
                
                % check added by sag
                if results(current_slice_number,left_right)~=0;
                    keyboard
                end
    
                results(current_slice_number,left_right)=fse.pixdim(1)*20*(poInd-ngInd+4000-2000*srlm)/10001; % slice # and left/right added by sag
                if show > 0
                    plot(20*(ngInd+1000+2000*srlm)/10001-10,ngVal+mean(yyShift),'x')
                    plot(20*(poInd+5000)/10001-10,poVal+mean(yyShift),'x')
                    title(['Left/Right: ' num2str(left_right)  '. Slice: ' num2str(current_slice_number)]); % added by sag 12dec16
                    hold off
                end

            end
        end
    end % for left_right, added by sag
    
    switch srlm
        case 0
            sp_out = results;
        case 1
            srlm_out = results;
    end
   
end

warning('on')
