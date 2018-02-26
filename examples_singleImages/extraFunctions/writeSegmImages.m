function writeSegmImages( config, fileSet, seg, segmentedCells)        

img = seg.image;


%% produce final segmentation mask

% small objects are excluded
% objects hiting the boundary are excluded


fSeeds = segmentedCells;

tsegmMask=zeros(size(img));
for nsc=1:numel( fSeeds )
    cellPixelInd = fSeeds(nsc).cellPixelListLindx;
    tsegmMask(cellPixelInd)=nsc;
end
L= tsegmMask;

%print
dpi=300;
dpi1 = 20;
fnameFa=[fileSet.resultsDir 'final_mask' num2str(fileSet.frameIdx) '.png'];
offScreenFig = figure('visible','off');
%figure;L=segmMask;
rgb = label2rgb(L,'jet',[0 0 0],'shuffle');
imshow(rgb,'Border', 'tight')
[H,W,~] = size(rgb);
set(offScreenFig, 'paperposition', [0 0 W/dpi1 H/dpi1]);
set(offScreenFig, 'papersize', [W/dpi1 H/dpi1]);
print(offScreenFig, '-dpng', sprintf('-r%d',dpi), fnameFa);
fprintf('Wrote %s\n', fnameFa);
close(offScreenFig);



% produce contour mask
segmMask=L;
se1=strel('square',3);linearIndAll=[];
for nsc = 1:max(segmMask(:))
    currObjMask = segmMask==nsc;
    currObjMaskDil = imdilate(currObjMask,se1);
    surPixMask=currObjMaskDil-currObjMask;    
    [Rv ,Cv] = find(surPixMask);    
    linearIndCur = sub2ind(size(img), Rv, Cv);
    % take the ones that are non-zero in the mask
    linearIndCurVals = segmMask(linearIndCur);
    linearIndCurValsAcc = linearIndCurVals>nsc; 
    linearIndAll= [linearIndAll;linearIndCur(linearIndCurValsAcc)];
end
linearIndCont = unique(linearIndAll);
segmMask(linearIndCont)=0;

offScreenFig = figure('visible','off');
fnameFb=[fileSet.resultsDir 'final_contour' num2str(fileSet.frameIdx) '.png'];
%figure;
imshow(img,[],'Border', 'tight')
hold on
% Display the initial contour on the original image in red.
contour(segmMask,[0  0],'r','LineWidth', 2)
set(offScreenFig, 'paperposition', [0 0 W/dpi1 H/dpi1]);
set(offScreenFig, 'papersize', [W/dpi1 H/dpi1]);
print(offScreenFig, '-dpng', sprintf('-r%d',dpi), fnameFb);
fprintf('Wrote %s\n', fnameFb);
close(offScreenFig);


% write control images
CellXPaperFigureWriter.writeSeedingControlImage( ...
    fileSet.oofImage, ...
    fileSet.seedingImageFile, ...
    seg.seeds, ...
    config...
    );


