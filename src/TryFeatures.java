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
import java.io.*;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/11/13
 * Time: 2:53 PM
 */
public class TryFeatures implements PlugInFilter, MouseListener {

    ImagePlus   inimg;
	ImagePlus 	profileImg;
    String      inimgPath;
    String      confFile;

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
		new ImagePlus("pattern", f.plotTemplate(new double[]{0, 0, 0})).show();

        inimg.show();

        inimg.getCanvas().addMouseListener(this);

        confFile = "current.conf";
        PrintWriter writer = null;
        try {
            writer = new PrintWriter(confFile);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        writer.print("");
        writer.close();

		profileImg = new ImagePlus();

    }

    public void mouseClicked(MouseEvent e)
    {

        ImageCanvas srcCanv = (ImageCanvas) e.getSource();
        int atX = 	srcCanv.offScreenX(e.getX());
		int atY = 	srcCanv.offScreenY(e.getY());

		FloatProcessor inip = (FloatProcessor) inimg.getProcessor().duplicate(); // make it possible to cast here  (FloatProcessor)

		//f.plotIntegResponse(atX, atY, inip);

		long t1, t2;
		t1 = System.currentTimeMillis();
        ImageProcessor ipPlot = f.getAngles(atX, atY, inip, true);
		t2 = System.currentTimeMillis();

		profileImg.setProcessor(ipPlot);
		profileImg.updateAndDraw();
		profileImg.show();

		IJ.log("\nextracted angles in "+((float)(t2-t1))+" msec.");

		if (f.ap!=null) {
			IJ.log("final angles:");
			for (int g = 0; g<f.ap.length; g++) {
				IJ.log("angle"+g+" > "+(f.ap[g]*(360f/((float)Math.PI*2)))+" degrees");
			}
		}

		t1 = System.currentTimeMillis();
		f.regionScores(atX, atY, inip, f.ap);
		t2 = System.currentTimeMillis();
		IJ.log("extracted features in "+((float)(t2-t1))+" msec.");

		if (f.ap!=null) {
			IJ.log("ap (sorted) angles:");
			for (int g = 0; g<f.ap.length; g++) {
				IJ.log("angle"+g+" > "+(f.ap[g]*(360f/((float)Math.PI*2)))+" degrees");
			}
		}

        double aA0 = (f.nA0>0)? (f.A0/f.nA0) : Double.NaN;
		IJ.log("A0 : "+f.A0+" ("+f.nA0+"), "+IJ.d2s(aA0, 1));

		IJ.log("---");

        double aA1 = (f.nA1>0)? (f.A1/f.nA1) : Double.NaN;
        IJ.log("A1 : "+f.A1+" ("+f.nA1+"),"+IJ.d2s(aA1, 1));
        double aA2 = (f.nA2>0)? (f.A2/f.nA2) : Double.NaN;
		IJ.log("A2 : "+f.A2+" ("+f.nA2+"),"+IJ.d2s(aA2, 1));
        double aA3 = (f.nA3>0)? (f.A3/f.nA3) : Double.NaN;
		IJ.log("A3 : "+f.A3+" ("+f.nA3+"),"+IJ.d2s(aA3, 1));

        IJ.log("---");

        double aB1 = (f.nB1>0)? (f.B1/f.nB1) : Double.NaN;
        IJ.log("B1 : "+f.B1+" ("+f.nB1+"),"+IJ.d2s((f.nB1>0)? (f.B1/f.nB1) : Double.NaN, 1));
        double aB2 = (f.nB2>0)? (f.B2/f.nB2) : Double.NaN;
		IJ.log("B2 : "+f.B2+" ("+f.nB2+"),"+IJ.d2s((f.nB2>0)? (f.B2/f.nB2) : Double.NaN, 1));
        double aB3 = (f.nB3>0)? (f.B3/f.nB3) : Double.NaN;
		IJ.log("B3 : "+f.B3+" ("+f.nB3+"),"+IJ.d2s((f.nB3>0)? (f.B3/f.nB3) : Double.NaN, 1));



		IJ.log("---");
        double ev0  = aA0-aB1;
        double ev1  = aA0-aB2;
        double ev2  = aA0-aB3;

        double ev3  = aA1-aB1;
        double ev4  = aA2-aB2;
        double ev5  = aA3-aB3;

//		double ev6  = aA2-aB1;
//        double ev7  = aA3-aB3; //(aA3>aB3)? ((aA3-aB3)/aA3) : 0;
//        double ev8  = aA3-aB2; //(aA3>aB2)? ((aA3-aB2)/aA3) : 0;

        String printEvidences = ""+IJ.d2s(ev0, 3)+" , "+IJ.d2s(ev1, 3)+" , "+IJ.d2s(ev2, 3)+" , "+
                IJ.d2s(ev3, 3)+" , "+IJ.d2s(ev4, 3)+" , "+
                IJ.d2s(ev5, 3);
//										+" , "+IJ.d2s(ev6, 1);
//		+" , "+
//                IJ.d2s(ev7, 1)+" , "+IJ.d2s(ev8, 1);

        IJ.log(printEvidences);

        // append line
        try {
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(confFile, true)));
            out.println(printEvidences);
            out.close();
        } catch (IOException e1) {
            //oh noes!
        }

		boolean showFeatures = false;
		if (showFeatures && f.ap!=null) {

			double[] ang_x = new double[9];
            double[] ang_y = new double[9];

            for (int i=0; i<9; i++) {
                ang_x[i] = i;
            }

			ang_y[0]   = ev0;
            ang_y[1]   = ev1;
            ang_y[2]   = ev2;
            ang_y[3]   = ev3;
            ang_y[4]   = ev4;
            ang_y[5]   = ev5;
//            ang_y[6]   = ev6;
//            ang_y[7]   = ev7;
//            ang_y[8]   = ev8;

			Plot p = new Plot("", "", "", ang_x, ang_y);
			p.draw();
			p.setColor(Color.RED);
			p.addPoints(ang_x, ang_y, Plot.BOX);
			p.show();

		}

        if (f.ap.length>=3) {

            ImagePlus templateFit = new ImagePlus("template", f.exportTemplate(f.ap));
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

		IJ.selectWindow("inimg");

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