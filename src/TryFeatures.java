import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/11/13
 * Time: 2:53 PM
 */
public class TryFeatures implements PlugInFilter, MouseListener {

    ImagePlus   inimg;
    String      inimgPath;
//	ImagePlus   inmask;
//    String      inmaskPath;

	Feat f;

    public void run(ImageProcessor imageProcessor)
	{

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
		//new ImagePlus("offsets",f.plotOffsets()).show();
		//new ImagePlus("pattern", f.exportTemplate(new double[]{0, Math.PI/2, Math.PI})).show();

        inimg.show();

        inimg.getCanvas().addMouseListener(this);

    }

    public void mouseClicked(MouseEvent e)
    {

        ImageCanvas srcCanv = (ImageCanvas) e.getSource();
        int atX = 	srcCanv.offScreenX(e.getX());
		int atY = 	srcCanv.offScreenY(e.getY());

		FloatProcessor inip = (FloatProcessor)inimg.getProcessor(); // make it possible to cast here

		//f.plotIntegResponse(atX, atY, inip);

		long t1, t2;
		t1 = System.currentTimeMillis();
        double[] angs = f.getAngles(atX, atY, inip, true);
		t2 = System.currentTimeMillis();

		IJ.log("\nextracted angles in "+((float)(t2-t1))+" msec.");

		if (angs!=null) {
			IJ.log("final angles:");
			for (int g = 0; g<angs.length; g++) {
				IJ.log("angle"+g+" > "+(angs[g]*(360f/((float)Math.PI*2)))+" degrees");
			}
		}

		t1 = System.currentTimeMillis();
		f.regionScores(atX, atY, inip, angs);
		t2 = System.currentTimeMillis();
		IJ.log("extracted features in "+((float)(t2-t1))+" msec.");

		if (f.ap!=null) {
			IJ.log("ap (sorted) angles:");
			for (int g = 0; g<f.ap.length; g++) {
				IJ.log("angle"+g+" > "+(f.ap[g]*(360f/((float)Math.PI*2)))+" degrees");
			}
		}

		IJ.log("A0 : "+f.A0+" ("+f.nA0+"),"+IJ.d2s((f.nA0>0)? (f.A0/f.nA0) : Double.NaN, 1));
		IJ.log("---");
		IJ.log("A1 : "+f.A1+" ("+f.nA1+"),"+IJ.d2s((f.nA1>0)? (f.A1/f.nA1) : Double.NaN, 1));
		IJ.log("A2 : "+f.A2+" ("+f.nA2+"),"+IJ.d2s((f.nA2>0)? (f.A2/f.nA2) : Double.NaN, 1));
		IJ.log("A3 : "+f.A3+" ("+f.nA3+"),"+IJ.d2s((f.nA3>0)? (f.A3/f.nA3) : Double.NaN, 1));
		IJ.log("---");
		IJ.log("B1 : "+f.B1+" ("+f.nB1+"),"+IJ.d2s((f.nB1>0)? (f.B1/f.nB1) : Double.NaN, 1));
		IJ.log("B2 : "+f.B2+" ("+f.nB2+"),"+IJ.d2s((f.nB2>0)? (f.B2/f.nB2) : Double.NaN, 1));
		IJ.log("B3 : "+f.B3+" ("+f.nB3+"),"+IJ.d2s((f.nB3>0)? (f.B3/f.nB3) : Double.NaN, 1));
		IJ.log("---");
		IJ.log(""+((f.A0/f.nA0)-(f.B1/f.nB1))/(f.A0/f.nA0));

		boolean showRegionAverages = false;
		if (showRegionAverages && angs!=null) {

			double[] ang_x = new double[7];
			ang_x[0]   = (f.nA0>3)? (f.A0/f.nA0) : (Double.MIN_VALUE);
			ang_x[1]   = (f.nA1>3)? (f.A1/f.nA1) : (Double.MIN_VALUE);
			ang_x[2]   = (f.nB1>3)? (f.B1/f.nB1) : (Double.MAX_VALUE);
			ang_x[3]   = (f.nA2>3)? (f.A2/f.nA2) : (Double.MIN_VALUE);
			ang_x[4]   = (f.nB2>3)? (f.B2/f.nB2) : (Double.MAX_VALUE);
			ang_x[5]   = (f.nA3>3)? (f.A3/f.nA3) : (Double.MIN_VALUE);
			ang_x[6]   = (f.nB3>3)? (f.B3/f.nB3) : (Double.MAX_VALUE);

			double[] ang_y = new double[7];
			ang_y[0] = 0;
			ang_y[1] = ang_y[2] = f.ap[0];
			ang_y[3] = ang_y[4] = f.ap[1];
			ang_y[5] = ang_y[6] = f.ap[2];

			Plot p = new Plot("", "", "", ang_y, ang_x);
			p.draw();
			p.setColor(Color.RED);
			p.addPoints(ang_y, ang_x, Plot.BOX);
			p.show();

		}

        if (angs.length>=3) {

            ImagePlus templateFit = new ImagePlus("template", f.exportTemplate(angs));
			templateFit.show();
            IJ.selectWindow("inimg");
            IJ.run("Add Image...", "image=template x="+(atX-templateFit.getWidth()/2)+" y="+(atY-templateFit.getHeight()/2)+" opacity=50");
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

	public int setup(String s, ImagePlus imagePlus) {
        if(imagePlus==null) {IJ.showMessage("needs image opened"); return DONE; }
		inimg = Tools.convertToFloatImage(imagePlus);
		inimg.setTitle("inimg");
		inimgPath = imagePlus.getOriginalFileInfo().directory+imagePlus.getOriginalFileInfo().fileName;
		return DOES_8G+DOES_32+NO_CHANGES;
	}

}