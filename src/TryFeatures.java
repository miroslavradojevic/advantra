import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.plugin.PlugIn;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
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
public class TryFeatures implements PlugIn, MouseListener {

    ImagePlus   inimg;
    String      inimgPath;
	ImagePlus   inmask;
    String      inmaskPath;

	Feat f;

    public void run(String s) {

        inimgPath = MyOpener.open("Open image file", false);
        inimg = new ImagePlus(inimgPath);
        if (inimg.getType()!=ImagePlus.GRAY32) {
            inimg = convertToFloatImage(inimg);
        }
        inimg.setTitle("inimg");

        inmaskPath = MyOpener.open("Open mask file", false);
        inmask = new ImagePlus(inmaskPath);
        inmask.setTitle("inmask");

        int     t       = 4;
        double  scale   = 2.0;
        double  D       = 10;
        GenericDialog gd = new GenericDialog("Fit Features");
        gd.addMessage("feature parameters");
        gd.addNumericField("neuron diameter", t, 0, 5, "pix");
        gd.addNumericField("check range radius", scale, 1, 5, "x diameter");
        gd.showDialog();
        if (gd.wasCanceled()) return;
		t 		=  	(int)gd.getNextNumber();
        scale   =   gd.getNextNumber();

		f= new Feat(t, scale);

//        inimg = convertToFloatImage(IJ.openImage());
//        inimg.setTitle("input_image");
//        ipFit = c.fit((FloatProcessor)inimg.getProcessor());

        inimg.show();
        inmask.show();
        IJ.selectWindow("inimg");
        IJ.run("Add Image...", "image=inmask x="+0+" y="+0+" opacity=30");
        //inmask.close();

        inimg.getCanvas().addMouseListener(this);

		IJ.log("filtering...");
		long t1 = System.currentTimeMillis();
        //new ImagePlus("", filterWithMask(inmask.getProcessor(), (FloatProcessor) inimg.getProcessor(), D)).show();
		long t2 = System.currentTimeMillis();
        IJ.log("done. "+((t2-t1)/1000f)+" sec.");

    }

    public FloatProcessor filterWithMask(ImageProcessor msk, FloatProcessor input, double D)
    {

        FloatProcessor ipOut = new FloatProcessor(input.getWidth(), input.getHeight());

        for (int x=0; x<input.getWidth(); x++) {
            for (int y=0; y<input.getHeight(); y++) {
                if (msk.getf(x, y)==255) {
                    double sc = f.bifurcationess(x, y, input, D);
                    ipOut.setf(x, y, (float)sc );
                }
            }
        }
        return ipOut;

    }

    public void mouseClicked(MouseEvent e)
    {

        ImageCanvas srcCanv = (ImageCanvas) e.getSource();
        int atX = 	srcCanv.offScreenX(e.getX());
		int atY = 	srcCanv.offScreenY(e.getY());

		FloatProcessor inip = (FloatProcessor)inimg.getProcessor(); // make it possible to cast here

        /*
		    visualization
		*/

		f.plotProfilesWithMSDetection(atX, atY,	inip);

        double[] angs = f.getAngles(atX, atY, inip, false);
        if (angs.length>=3) {

            ImagePlus templateFit = new ImagePlus("template", f.exportTemplate(angs));
			templateFit.show();
            IJ.selectWindow("inimg");
            IJ.run("Add Image...", "image=template x="+(atX-templateFit.getWidth()/2)+" y="+(atY-templateFit.getHeight()/2)+" opacity=20");
			templateFit.close();

            /*
            ImagePlus showC;
            showC = new ImagePlus("feature.template", f.exportTemplate(angs));
            showC.show();
            for (int q=0; q<5; q++) showC.getCanvas().zoomIn(0, 0);
            */

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

//	public int setup(String s, ImagePlus imagePlus) {
//        if(imagePlus==null) {IJ.showMessage("needs image to work!"); return DONE; }
//		inimg = convertToFloatImage(imagePlus);
//		inimg.setTitle("input_image");
//		return DOES_8G+DOES_32+NO_CHANGES;
//	}

}