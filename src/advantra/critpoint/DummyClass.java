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
 * To change this template use File | Settings | File Templates.
 */
public class DummyClass implements PlugIn, MouseListener {

    Conf3 c3;
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
        d =  (int)gd.getNextNumber();
        Xon =  (float)gd.getNextNumber();
        Xoff =  (float)gd.getNextNumber();

        System.out.println("loaded -> "+d+" , "+Xon+" and "+Xoff);

        c3 = new Conf3(d);
		new ImagePlus("all", c3.plotAll()).show();
		//new ImagePlus("one", c3.plotOne(37)).show();


        ImageStack is = new ImageStack(img.getWidth(), img.getHeight());
        for (int t = 0; t < c3.regionIdxMap.size(); t++) {

            System.out.println("calculating feature " + t + " on whole image...");

            ImageProcessor ip2 = new FloatProcessor(img.getWidth(), img.getHeight());

            // compare with calculation per point
            for (int x = 0; x < img.getWidth(); x++) {
                for (int y =0; y < img.getHeight(); y++) {

                    float val = (float) c3.scoreAtPos(x, y, t, img.getProcessor(), Xon, Xoff);
                    ip2.setf(x, y, val);

                }
            }

            is.addSlice(ip2);

        }

		ImagePlus      im2 = new ImagePlus("score.all", is);
		im2.show();

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

    @Override
    public void mouseClicked(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.


        ImageCanvas srcCanv = (ImageCanvas) e.getSource();

        String source  =  srcCanv.getName();

		int atX = 	srcCanv.offScreenX(e.getX());
		int atY = 	srcCanv.offScreenY(e.getY());

        System.out.println("X : "+atX+"    ,     Y : "+atY);
        float val = (float) c3.scoreAtPos(atX, atY, 37, img.getProcessor(), 80, 50);

    }

    @Override
    public void mousePressed(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void mouseReleased(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void mouseEntered(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void mouseExited(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }
}
