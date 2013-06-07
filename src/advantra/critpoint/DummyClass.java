package advantra.critpoint;

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
 * Date: 5/31/13
 * Time: 10:12 AM
 */
public class DummyClass implements PlugIn, MouseListener {

    Conf3 c3;
    CircularConfiguration3 ccf3;
    ImagePlus img;

	public void run(String s) {

		img = convertToFloatImage(IJ.openImage());
		img.show();

        ImageCanvas imgCanv = img.getCanvas();
        imgCanv.addMouseListener(this);

		int d = 10;//img.getHeight();
        float Xon = 80;
        float Xoff = 40;

        GenericDialog gd = new GenericDialog("feature extraction parameters");
        gd.addNumericField("patch radius", d, 0, 5, "pixels");
		gd.addNumericField("XON",  Xon, 0, 5, " defines exponential curve for ON");
        gd.addNumericField("XOFF", Xoff, 0, 5, " defines exponential curve for OFF");
        gd.showDialog();
        if (gd.wasCanceled()) return;
        d       =  (int)gd.getNextNumber();
        Xon     =  (float)gd.getNextNumber();
        Xoff    =  (float)gd.getNextNumber();

/*        c3 = new Conf3(d);
		new ImagePlus("all", c3.plotAll()).show();

        ImageStack is = new ImageStack(img.getWidth(), img.getHeight());
        for (int t = 0; t < c3.regionIdxMap.size(); t++) {

            ImageProcessor ip2 = new FloatProcessor(img.getWidth(), img.getHeight());

            // compare with calculation per point
            for (int x = d; x < img.getWidth()-d; x++) {
                for (int y = d; y < img.getHeight()-d; y++) {

                    float val = 1;//(float) c3.scoreAtPos(x, y, t, img.getProcessor(), Xon, Xoff);
                    ip2.setf(x, y, val);

                }
            }

            is.addSlice(ip2);

        }

		ImagePlus      im2 = new ImagePlus("score.all.multiplication", is);
		im2.show();*/


        /*
        ccf3
         */



        ccf3 = new CircularConfiguration3(d);
        ImagePlus feats3 = new ImagePlus("Y.feat."+d, ccf3.plotKernels());
        feats3.show();
        feats3.getCanvas().zoomIn(0, 0);

        for (int t = 0; t < 2; t++) { // ccf3.kernels.size()

//            ImageStack isExperimental = new ImageStack(img.getWidth(), img.getHeight());
//            ImageProcessor ipExperimental = new FloatProcessor(img.getWidth(), img.getHeight());
//            isExperimental.addSlice(ipExperimental);


              System.out.println("scoring for "+t+" feature...");
              ImageStack isExperimental = ccf3.score_Experimental(t, img.getProcessor());

//            // compare with calculation per point
//            for (int x = d; x < img.getWidth()-d; x++) {
//                for (int y = d; y < img.getHeight()-d; y++) {
//
//                    float val = (float) c3.scoreAtPos(x, y, t, img.getProcessor(), Xon, Xoff);
//                    ip2.setf(x, y, val);
//
//                }
//            }

            ImagePlus      imExperimental = new ImagePlus("score.ap_an."+t+"feature", isExperimental);
            imExperimental.show();

        }



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

    public void mouseClicked(MouseEvent e) {

        ImageCanvas srcCanv = (ImageCanvas) e.getSource();

        String source  =  srcCanv.getName();

		int atX = 	srcCanv.offScreenX(e.getX());
		int atY = 	srcCanv.offScreenY(e.getY());

        System.out.println("X : "+atX+"    ,     Y : "+atY);
        float val = (float) c3.scoreAtPos(atX, atY, 37, img.getProcessor(), 80, 50);

    }

    public void mousePressed(MouseEvent e) {
    }

    public void mouseReleased(MouseEvent e) {
    }

    public void mouseEntered(MouseEvent e) {
    }

    public void mouseExited(MouseEvent e) {
    }
}
