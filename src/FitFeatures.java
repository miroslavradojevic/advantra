import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/11/13
 * Time: 2:53 PM
 */
public class FitFeatures implements PlugIn, MouseListener {

    ImageProcessor ipFit;
    ImagePlus inimg;
	ImagePlus profileImage;

	Feat f;

    public void run(String s) {

		//default
        int t = 5;
        double scale = 2.0;

        GenericDialog gd = new GenericDialog("Fit Features");
        gd.addMessage("feature parameters");
        gd.addNumericField("neuron diameter", t, 0, 5, "pix");
        gd.addNumericField("check range radius", scale, 0, 5, "x diameter");
        gd.showDialog();
        if (gd.wasCanceled()) return;
		t 		=  	(int)gd.getNextNumber();
        scale   =   gd.getNextNumber();

		f= new Feat(t, scale);

        ImagePlus showC;
        showC = new ImagePlus("", f.plotOffsets());
        /*showC.show();
        for (int q=0; q<5; q++) showC.getCanvas().zoomIn(0, 0);*/

        //IJ.showMessage("Open image to fit the configurations on");
        inimg = convertToFloatImage(IJ.openImage());
        inimg.setTitle("input_image");

//        IJ.showMessage("start calculating best configurations of the feature...");
//        IJ.log("wait, calculating...");
//        ipFit = c.fit((FloatProcessor)inimg.getProcessor());
//        IJ.log("finished.");
//        IJ.showMessage("found best fits");

        inimg.show();
        inimg.getCanvas().addMouseListener(this);
        //IJ.showMessage("click on the location to see the configuration fit there");

    }

    public void mouseClicked(MouseEvent e) {

        ImageCanvas srcCanv = (ImageCanvas) e.getSource();
        int atX = 	srcCanv.offScreenX(e.getX());
		int atY = 	srcCanv.offScreenY(e.getY());

		FloatProcessor inip = (FloatProcessor)inimg.getProcessor();
		int Npoints = 200;
		int maxIterBlurring = 15;
		int maxIterNonBlurring = 200;
		double eps = 0.0001;
		double minD = 0.2;
		int M = 1;


        double[] angs = f.getAngles(atX, atY, inip, Npoints, maxIterNonBlurring, eps, M, minD);
		double[] sc = f.scores(atX, atY, inip, Npoints, maxIterNonBlurring, eps, M, minD);

		/*
			visualization
		 */

		ImageStack is = f.plotProfilesWithMSDetection(atX, atY,	inip, Npoints, maxIterBlurring,	maxIterNonBlurring,	eps, minD, M);

		if (profileImage==null) {     // create it if it did not exist
			profileImage = new ImagePlus("profile", is);
			profileImage.show();
		}
		else {
			profileImage.setStack(is);
			profileImage.updateImage();
			profileImage.draw();

		}

        if (angs!=null) {

			ImagePlus templateFit = new ImagePlus("template", f.exportTemplate(angs));
			templateFit.show();
            IJ.selectWindow("input_image");
            IJ.run("Add Image...", "image=template x="+(atX-templateFit.getWidth()/2)+" y="+(atY-templateFit.getHeight()/2)+" opacity=50");
			templateFit.close();

		}

    }

    public void mousePressed(MouseEvent e) {
    }

    public void mouseReleased(MouseEvent e) {
    }

    public void mouseEntered(MouseEvent e) {
    }

    public void mouseExited(MouseEvent e) {
    }

    private ImagePlus convertToFloatImage(
            ImagePlus inim
    )
    {

        int W = inim.getWidth();
        int H = inim.getHeight();

        ImageProcessor ip = new FloatProcessor(W, H);
        for (int i = 0; i < H*W; i++ ) {
            ip.setf(i, inim.getProcessor().getPixelValue(i%W, i/W));
        }

        ImagePlus outim = new ImagePlus(inim.getTitle()+"ToFloat", ip);
        return outim;

    }

}