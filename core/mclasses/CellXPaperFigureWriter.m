classdef CellXPaperFigureWriter < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        
        function writeInitialInputImage(srcImgFile, dstImgFile, seeds, config)
            img = CellXImageIO.loadToGrayScaleImage(srcImgFile, config.cropRegionBoundary);
            img = CellXImageIO.normalize( img );
            %img = CellXImgGen.generateHoughTransformImage(seeds, img);
            imwrite(img, dstImgFile, 'png');
            fprintf('Wrote %s\n', dstImgFile);
        end
        
        
        function writeHoughTransformControlImage(srcImgFile, dstImgFile, seeds, config)
            img = CellXImageIO.loadToGrayScaleImage(srcImgFile, config.cropRegionBoundary);
            img = CellXImageIO.normalize( img );
            img = CellXImgGen.generateHoughTransformImage(seeds, img);
            imwrite(img, dstImgFile, 'png');
            fprintf('Wrote %s\n', dstImgFile);
        end
        
        function writeSegmentationControlImage(srcImgFile, dstImgFile, seeds, config)
            img = CellXImageIO.loadToGrayScaleImage(srcImgFile, config.cropRegionBoundary);
            img = CellXImageIO.normalize( img );
            img = CellXImgGen.generateControlImage(seeds, img);
            imwrite(img, dstImgFile, 'png');
            fprintf('Wrote %s\n', dstImgFile);
        end
        
        function writeSegmentationControlImageWithIndices(srcImgFile, dstImgFile, seeds, config)
            img = CellXImageIO.loadToGrayScaleImage(srcImgFile, config.cropRegionBoundary);
            img = CellXImageIO.normalize( img );
            CellXImgGen.createControlImageWithCellIndices(seeds, img,  dstImgFile);
        end
        
        function writeSeedingControlImage(srcImgFile, dstImgFile, seeds, config)
            img = CellXImageIO.loadToGrayScaleImage(srcImgFile, config.cropRegionBoundary);
            img = CellXImageIO.normalize( img );
            img = CellXImgGen.generateHoughTransformImage(seeds, img);
            imwrite(img, dstImgFile, 'png');
            fprintf('Wrote %s\n', dstImgFile);
        end
        
       function writeSeedingControlImageComparison(srcImgFile, dstImgFile, seeds, config)
            img = CellXImageIO.loadToGrayScaleImage(srcImgFile, config.cropRegionBoundary);
            img = CellXImageIO.normalize( img );
            [rgb imgNuc] = CellXImgGen.generateHoughTransformImageComparison(seeds, img);
            [fpathstr, fname] = fileparts(dstImgFile);
            rgbfilename = [fpathstr filesep 'RGB' fname];
            nucfilename = [fpathstr filesep 'Mask' fname];
            imwrite(rgb, [rgbfilename '.tiff'], 'tiff');
            imwrite(double(imgNuc), [nucfilename '.tiff'], 'tiff');
            fprintf('Wrote %s\n', dstImgFile);
        end
        
        function writeSegmentationControlImageWithIndicesB(srcImgFile, dstImgFile, seeds, config)
            img = CellXImageIO.loadToGrayScaleImage(srcImgFile, config.cropRegionBoundary);
            img = CellXImageIO.normalize( img );
            CellXImgGen.createControlImageWithCellIndicesB(seeds, img,  dstImgFile);
        end
        
        
        function writeFluoControlImage(srcImgFile, dstImgFile, seeds, config)
            img = CellXImageIO.loadToGrayScaleImage(srcImgFile, config.cropRegionBoundary);
            img = CellXImageIO.normalize( img );
            img = CellXImgGen.generateIntensityExtractionControlImage( ...
                seeds, ...
                img ...
                );
            imwrite(img, dstImgFile, 'png');
            fprintf('Wrote %s\n', dstImgFile);
        end
        
        function writeFluoControlBoundaryImage(srcImgFile, dstImgFile, seeds, config)
            img = CellXImageIO.loadToGrayScaleImage(srcImgFile, config.cropRegionBoundary);
            img = CellXImageIO.normalize( img );
            img = CellXImgGen.generateIntensityExtractionControlBoundaryImage( ...
                seeds, ...
                img ...
                );
            imwrite(img, dstImgFile, 'png');
            fprintf('Wrote %s\n', dstImgFile);
        end
        
        function writeGradientControlBoundaryImage(srcImgFile, dstImgFile, seeds, config)
            img = CellXImageIO.loadToGrayScaleImage(srcImgFile, config.cropRegionBoundary);           
            [grdX grdY] = gradient(img);
            timageGradientMagn = sqrt(grdX.^2 + grdY.^2);
            img = CellXImageIO.normalize( timageGradientMagn );
            img = CellXImgGen.generateIntensityExtractionControlBoundaryImage( ...
                seeds, ...
                img ...
                );
            imwrite(img, dstImgFile, 'png');
            fprintf('Wrote %s\n', dstImgFile);
        end
        
        
        
        function writeLineageImage(dstImgFile, trackInfoMatrix)
            CellXImgGen.generateLineageTree(trackInfoMatrix, dstImgFile);
            fprintf('Wrote %s\n', dstImgFile);
        end 
        
        
    end
    
end

