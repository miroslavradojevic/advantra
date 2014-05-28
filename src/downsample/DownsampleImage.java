package downsample;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.awt.*;
import java.awt.event.TextEvent;
import java.awt.event.TextListener;
import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/15/13
 * Time: 10:59 AM
 */
public class DownsampleImage implements PlugInFilter {

	/*
	plugin that downsamples the image with given parameters
	- new width in pixels (height is scaled respectively)
	- sourceSigma
	- targetSigma
	 */

    //static int MIN_W=1, MIN_H=1;
	//static int NR_PIX_PER_DIAM = 4;

    ImagePlus imp;

    int     width, height;
    int     newWidth, newHeight;
    float   sourceSigma, targetSigma;
    //boolean keepSource;
    GenericDialog gd;

    TextListener tl = new TextListener() {
        public void textValueChanged(TextEvent e) {
            TextField tf = (TextField)e.getSource();
            if (!tf.getText().equals("")){
                // new value was filled
                if ( tf.getName().equals("widthPrompt") ) {
                    newWidth    = Integer.valueOf(tf.getText());
                    // set height accordingly
                    if (newWidth<width && newWidth>0) {
                        newHeight   = Math.round( newWidth * height / width );
                    }
                    else {

                        // reset (irregular values at the input)
                        newHeight = height;
                        newWidth = width;

                    }

                    gd.getTextArea1().setText(String.valueOf(newWidth));
                    gd.getTextArea2().setText(String.valueOf(newHeight));

                }
            }
            else {
                // it was empty string - make both empty
                newWidth = width;
                newHeight = height;

                gd.getTextArea1().setText("");
                gd.getTextArea2().setText("");

            }
        }
    };

    public int setup(String s, ImagePlus imagePlus) {

        this.imp = imagePlus;

        width = imagePlus.getWidth();
        height = imagePlus.getHeight();

		newWidth = width;
		newHeight = height;

        sourceSigma = .5f;
        targetSigma = .5f;

        gd = new GenericDialog("DOWNSAMPLE IMAGE");
		gd.addMessage(width+" x "+height);
		gd.addNumericField(	"width  :",  newWidth,  0, 5, " pix.");
		gd.addTextAreas("newWidth", "newHeight", 1, 10);

        // add listener to width field
        TextField promptWidth = (TextField)gd.getNumericFields().get(0);
        promptWidth.setName("widthPrompt");
        promptWidth.addTextListener(tl);

        gd.addNumericField("sourceSigma",        sourceSigma, 1, 10, " pix.");
        gd.addNumericField("targetSigma",        targetSigma, 1, 10, " pix.");

        gd.showDialog();

        if (gd.wasCanceled()) return DONE;

        //width       = (int) gd.getNextNumber();
        //height      = (int) gd.getNextNumber();

        return DOES_8G+NO_CHANGES;
    }

    public void run(ImageProcessor imageProcessor) {

        System.out.println("downsampling to: " + newWidth + " x " + newHeight);

        System.out.println(
                "CALIBRATION:\n" +
                    "pixel = \n" + imp.getCalibration().pixelWidth + " x " + imp.getCalibration().pixelHeight +"" +
                        "\n["+ imp.getCalibration().getUnit()+"]");
        System.out.println("sourceSigma: " + sourceSigma + "\ntargetSigma: " + targetSigma);

        ImagePlus imDown = downsample(imp, newWidth, sourceSigma, targetSigma);
        imDown.show();

        System.out.println(
                "CALIBRATION:\n" +
                     "pixel = \n" + imDown.getCalibration().pixelWidth + " x " + imDown.getCalibration().pixelHeight +"" +
                        "\n["+ imp.getCalibration().getUnit()+"]");


    }

    public ImagePlus downsample(ImagePlus inimg, int newWidth, float sourceSigma, float targetSigma){

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
