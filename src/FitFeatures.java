import ij.IJ;
import ij.ImagePlus;
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

    Conf c;

	Feat f;

    public void run(String s) {

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

        //c = new Conf(t, scale);

		if(true) return;

//        System.out.println("\n\n SUMMARY: \n\n");
//        System.out.println(""+c.regionIdxMap.size()+" , "+c.regionSize.size()+" index maps");

//        for (int m = 0; m<c.regionIdxMap.size(); m++) {
//
//            System.out.println("\nReg. "+m+" : rotation "+(m/Conf.nRot));
//            System.out.println(""+c.names.get(m/Conf.nRot)+" ");
//            for (int i = 0; i<c.regionSize.get(m).length; i++) {
//                System.out.print("["+i+" -> "+c.regionSize.get(m)[i]+"],");
//            }
//
//        }

//        System.out.println(""+c.angles.size()+" angle");
//        System.out.println(""+c.names.size()+" names");
//        System.out.println("\n\n --- \n\n");

        ImagePlus showC;

        showC = new ImagePlus("templates(diam="+c.diam+",r="+c.r+")", c.plotTemplatesAll());

        showC.show();
        showC.getCanvas().zoomIn(0, 0);
        showC.getCanvas().zoomIn(0, 0);
        showC.getCanvas().zoomIn(0, 0);
        showC.getCanvas().zoomIn(0, 0);


        showC = new ImagePlus("templates(diam="+c.diam+",r="+c.r+")", c.plotTemplates());

        showC.show();
        showC.getCanvas().zoomIn(0, 0);
        showC.getCanvas().zoomIn(0, 0);
        showC.getCanvas().zoomIn(0, 0);
        showC.getCanvas().zoomIn(0, 0);


















        //if(true) return;

        showC = new ImagePlus("kernels(diam="+c.diam+",r="+c.r+")", c.plotKernels());

//        showC.show();
//        showC.getCanvas().zoomIn(0, 0);
//        showC.getCanvas().zoomIn(0, 0);
//        showC.getCanvas().zoomIn(0, 0);
//        showC.getCanvas().zoomIn(0, 0);

        // fit
        IJ.showMessage("Open image to fit the configurations on");
        ImagePlus inimg = convertToFloatImage(IJ.openImage());
        inimg.setTitle("input_image");

        IJ.showMessage("start calculating best configurations of the feature...");
        IJ.log("wait, calculating...");
        ipFit = c.fit((FloatProcessor)inimg.getProcessor());
        IJ.log("finished.");
        IJ.showMessage("found best fits");

        inimg.show();
        inimg.getCanvas().addMouseListener(this);
        IJ.showMessage("click on the location to see the configuration fit there");

    }

    public void mouseClicked(MouseEvent e) {

        ImageCanvas srcCanv = (ImageCanvas) e.getSource();
        int atX = 	srcCanv.offScreenX(e.getX());
		int atY = 	srcCanv.offScreenY(e.getY());

        System.out.println("X: "+atX+" , Y: "+atY);

        int configurationIdx = (int)Math.round(ipFit.getf(atX, atY));
        System.out.println("configuration index to plot here: "+configurationIdx);
        // show one that fits there

        ImagePlus bestFitImage = new ImagePlus("a1", c.plotKernel(configurationIdx));
        bestFitImage.show();
        //System.out.println(bestFitImage.getShortTitle());
        IJ.selectWindow("input_image");
        IJ.run("Add Image...", "image=a1 x="+(atX-bestFitImage.getWidth()/2)+" y="+(atY-bestFitImage.getHeight()/2)+" opacity=50");
        bestFitImage.close();

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