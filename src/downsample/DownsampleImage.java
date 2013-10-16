package downsample;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
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
public class DownsampleImage implements PlugInFilter {      // , TextListener

	/*
	plugin that downsamples the image with given parameters
	- new width in pixels (height is scaled respectively)
	- sourceSigma
	- targetSigma
	-
	 */

    static int MIN_W=10, MIN_H=10;
	static int NR_PIX_PER_DIAM = 4;

    int     width, height;
    int     newWidth, newHeight;
    float   sourceSigma, targetSigma;
    boolean keepSource;
    GenericDialog gd;

    TextListener tl = new TextListener() {
        public void textValueChanged(TextEvent e) {

            System.out.println("changed!");


            if(true){

            TextField tf = (TextField)e.getSource();

            if (!tf.getText().equals("")){

                // new value was filled
                if ( tf.getName().equals("widthPrompt") ) {


                    System.out.println("width changed!");

                    newWidth    = Integer.valueOf(tf.getText());
                    TextField vc;

                    // set height accordingly
                    if (newWidth<width && newWidth>MIN_W) {
                        newHeight   = Math.round( newWidth * height / width );
//                        vc = (TextField)gd.getNumericFields().get(1);
//                        vc.setText(String.valueOf(newHeight));
//                        gd.getNumericFields().set(1, vc);
                    }
                    else {
                        // reset, fake values at the input
                        newHeight = height;
                        newWidth = width;

//                        vc = (TextField)gd.getNumericFields().get(0);
//                        vc.setText(String.valueOf(width));
//                        gd.getNumericFields().set(0, vc);

//                        vc = (TextField)gd.getNumericFields().get(1);
//                        vc.setText(String.valueOf(height));
//                        gd.getNumericFields().set(1, vc);

                    }

                }

//                if ( tf.getName().equals("heightPrompt") ) {
//
//                    System.out.println("height changed!");
//
//                    newHeight    = Integer.valueOf(tf.getText());
//                    TextField vc;
//
//                    // set width accordingly
//                    if (newHeight<height && newHeight>MIN_H) {
//                        newWidth   = Math.round( newHeight * width / height );
//                        vc = (TextField)gd.getNumericFields().get(1);
//                        vc.setText(String.valueOf(newWidth));
//                        gd.getNumericFields().set(0, vc);
//                    }
//                    else {
//                        // reset, fake values at the input
//                        newHeight = height;
//                        newWidth = width;
//
//                        vc = (TextField)gd.getNumericFields().get(0);
//                        vc.setText(String.valueOf(width));
//                        gd.getNumericFields().set(0, vc);
//
//                        vc = (TextField)gd.getNumericFields().get(1);
//                        vc.setText(String.valueOf(height));
//                        gd.getNumericFields().set(1, vc);
//
//                    }
//
//                }

            }
            else {
                // it was empty string - make both empty

                System.out.println("empty!");
                newWidth = width;
                newHeight = height;

//                TextField vc;
//                vc = (TextField)gd.getNumericFields().get(0);
//                vc.setText("");
//                gd.getNumericFields().set(0, vc);
//                vc = (TextField)gd.getNumericFields().get(1);
//                vc.setText("");
//                gd.getNumericFields().set(1, vc);

            }

            }






        }
    };

    public int setup(String s, ImagePlus imagePlus) {

        width = imagePlus.getWidth();
        height = imagePlus.getHeight();

		newWidth = width;
		newHeight = height;

        sourceSigma = .5f;
        targetSigma = .5f;

		System.out.println("pixel = \n" + imagePlus.getCalibration().pixelWidth + " x " + imagePlus.getCalibration().pixelHeight +"  ["+ imagePlus.getCalibration().getUnit()+"]");

        gd = new GenericDialog("DOWNSAMPLE IMAGE");
		gd.addMessage("W x H = "+width+" x "+height);
		gd.addMessage("--------");
		gd.addNumericField(	"width  :",  newWidth,  0, 5, " pix.");
        gd.addMessage( 		"height :" + newHeight);
		gd.addTextAreas("t1", "t2", 1, 1);

        // add listener to width field
        TextField promptWidth = (TextField)gd.getNumericFields().get(0);
        promptWidth.setName("widthPrompt");
        promptWidth.addTextListener(tl);

        //gd.addNumericField("height",        height, 0, 5, " pix.");
        // add listener to this field
        //TextField promptHeight = (TextField)gd.getNumericFields().get(1);
        //promptHeight.setName("heightPrompt");
        //promptHeight.addTextListener(tl);

        gd.showDialog();

        //prompt.setBackground(Color.BLUE);

        if (gd.wasCanceled()) return DONE;

        //width       = (int) gd.getNextNumber();
        //height      = (int) gd.getNextNumber();

        System.out.println("downsampling to: " + newWidth + " x " + newHeight);

        return DOES_8G+NO_CHANGES;
    }

    public void run(ImageProcessor imageProcessor) {




    }

//    public void textValueChanged(TextEvent e) {
//        IJ.log("text changed");
//    }
}
