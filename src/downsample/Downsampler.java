package downsample;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.Roi;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/15/13
 * Time: 10:57 AM
 */
public class Downsampler {

    public Downsampler(){}

    public ImagePlus run(ImagePlus inimg, int newWidth, float sourceSigma, float targetSigma){

        if (newWidth<=inimg.getWidth() && newWidth>0){

            // scale height
            int newHeight   = Math.round( newWidth * inimg.getHeight() / inimg.getWidth() );

            // make duplicate
            ImagePlus outimg = inimg.duplicate();
            outimg.setTitle("DOWNSAMPLED_"+inimg.getTitle()+"_"+newWidth+"x"+newHeight);

            float s = targetSigma * inimg.getWidth() / newWidth;
            IJ.log("sigma used = " + Math.sqrt( s * s - sourceSigma * sourceSigma ));
            IJ.run(outimg, "Gaussian Blur...", "sigma=" + Math.sqrt( s * s - sourceSigma * sourceSigma ) + " stack" );
            IJ.run(outimg, "Scale...", "x=- y=- width=" + newWidth + " height=" + newHeight + " process title=- interpolation=None" );

            float extraX = (inimg.getWidth() % 2 == 0) ? 0 : 1;
            float extraY = (inimg.getHeight() % 2 == 0) ? 0 : 1;
            float initialX = (newWidth % 2 == 0) ? (inimg.getWidth() / 2 - newWidth/2 + extraX) : (inimg.getWidth() / 2 - newWidth/2 +1 -extraX);
            float initialY = (newHeight % 2 == 0) ? (inimg.getHeight() / 2 - newHeight/2 + extraY) : (inimg.getHeight() / 2 - newHeight/2 +1 -extraY);

            outimg.setRoi(new Roi(initialX, initialY, newWidth, newHeight));
            IJ.run(outimg, "Crop", "");
//                IJ.makeRectangle(initialX, initialY, width, height);
            IJ.run(outimg, "Canvas Size...", "width=" + newWidth + " height=" + newHeight + " position=Center" );

            outimg.getCalibration().pixelWidth      *= ((double) outimg.getWidth()     / newWidth);
            outimg.getCalibration().pixelHeight     *= ((double) outimg.getHeight()    / newHeight);

            return outimg;

        }
        else{
            IJ.log("parameter set for upsampling");
            return null;
        }

    }

}
