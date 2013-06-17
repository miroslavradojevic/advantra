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
    ImagePlus inimg;

//    Conf c;

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

//        c = new Conf(t, scale);

//        for (int m = 0; m<c.regionIdxMap.size(); m++) {
//            System.out.println("\nReg. "+m+" : rotation "+(m/Conf.nRot));
//            System.out.println(""+c.names.get(m/Conf.nRot)+" ");
//            for (int i = 0; i<c.regionSize.get(m).length; i++) {
//                System.out.print("["+i+" -> "+c.regionSize.get(m)[i]+"],");
//            }
//        }

        ImagePlus showC;

//        showC = new ImagePlus("templates(diam="+c.diam+",r="+c.r+")", c.plotTemplatesAll());
//
//        showC.show();
//        showC.getCanvas().zoomIn(0, 0);
//        showC.getCanvas().zoomIn(0, 0);
//        showC.getCanvas().zoomIn(0, 0);
//        showC.getCanvas().zoomIn(0, 0);

//        showC = new ImagePlus("templates(diam="+c.diam+",r="+c.r+")", c.plotTemplates());
//
//        showC.show();
//        showC.getCanvas().zoomIn(0, 0);
//        showC.getCanvas().zoomIn(0, 0);
//        showC.getCanvas().zoomIn(0, 0);
//        showC.getCanvas().zoomIn(0, 0);

        showC = new ImagePlus("", f.showOffsets());

        showC.show();
        showC.getCanvas().zoomIn(0, 0);
        showC.getCanvas().zoomIn(0, 0);
        showC.getCanvas().zoomIn(0, 0);
        showC.getCanvas().zoomIn(0, 0);

        IJ.showMessage("Open image to fit the configurations on");
        inimg = convertToFloatImage(IJ.openImage());
        inimg.setTitle("input_image");

//        IJ.showMessage("start calculating best configurations of the feature...");
//        IJ.log("wait, calculating...");
//        ipFit = c.fit((FloatProcessor)inimg.getProcessor());
//        IJ.log("finished.");
//        IJ.showMessage("found best fits");

        inimg.show();
        inimg.getCanvas().addMouseListener(this);
        IJ.showMessage("click on the location to see the configuration fit there");

    }

    public void mouseClicked(MouseEvent e) {

        ImageCanvas srcCanv = (ImageCanvas) e.getSource();
        int atX = 	srcCanv.offScreenX(e.getX());
		int atY = 	srcCanv.offScreenY(e.getY());

        //int configurationIdx = (int)Math.round(ipFit.getf(atX, atY));
        double[] angs = f.get3Angles(atX, atY, (FloatProcessor)inimg.getProcessor(), 200, 200, 0.0001, 1, 0.5);

        angs = new double[]{0.0, 1.0, 2.0};
        //angs = new double[]{Math.PI, Math.PI+Math.PI/2, Math.PI*2};

        if (angs!=null) {

            for (int g = 0; g < angs.length; g++) {
                IJ.log("arg: "+angs[g]);
            }

            int[][] regionMap = new int[1][];
            int[][] regionSize = new int[1][];
            float[][] kernel = new float[1][];

            ImageProcessor templ = f.template(angs, regionMap, regionSize, kernel);
            ImagePlus bestFitImage = new ImagePlus("a1", templ);
            bestFitImage.show();
            IJ.selectWindow("input_image");
            IJ.run("Add Image...", "image=a1 x="+(atX-bestFitImage.getWidth()/2)+" y="+(atY-bestFitImage.getHeight()/2)+" opacity=50");
            bestFitImage.close();
            new ImagePlus("", f.plotSums(atX, atY, (FloatProcessor) inimg.getProcessor())).show();
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