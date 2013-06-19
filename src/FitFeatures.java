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
        int t = 4;
        double scale = 1.5;

        GenericDialog gd = new GenericDialog("Fit Features");
        gd.addMessage("feature parameters");
        gd.addNumericField("neuron diameter", t, 0, 5, "pix");
        gd.addNumericField("check range radius", scale, 1, 5, "x diameter");
        gd.showDialog();
        if (gd.wasCanceled()) return;
		t 		=  	(int)gd.getNextNumber();
        scale   =   gd.getNextNumber();

		f= new Feat(t, scale);

        /*
        ImagePlus showC;
        showC = new ImagePlus("", f.plotOffsets());
        showC.show();
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

        new ImagePlus("scores", f.scores((FloatProcessor)inimg.getProcessor())).show();

    }

    public void mouseClicked(MouseEvent e)
    {

        ImageCanvas srcCanv = (ImageCanvas) e.getSource();
        int atX = 	srcCanv.offScreenX(e.getX());
		int atY = 	srcCanv.offScreenY(e.getY());

		FloatProcessor inip = (FloatProcessor)inimg.getProcessor();

		double[] sc = f.scores(atX, atY, inip);    // TODO; add partial results
        if (sc!=null) IJ.log("\nscore : "+sc[0]+"\n");

		/*
			visualization
		 */
        int maxIterBlurring = 15;
		f.plotProfilesWithMSDetection(atX, atY,	inip);

        double[] angs = f.getAngles(atX, atY, inip, false);
        if (angs.length>=3) {
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