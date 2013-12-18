package detection2d;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/12/13
 * Time: 11:45 AM
 * usage of the class Masker2D, plugin
 */
public class Masker2DDemo implements PlugInFilter {

    float[][] 	inimg_xy; // store input image as an array

    /*
    parameters
     */
    float       nhoodRadius, iDiff;
    int         CPU_NR;

    public int setup(String s, ImagePlus imagePlus) {

		if(imagePlus==null) return DONE;

		inimg_xy = new float[imagePlus.getWidth()][imagePlus.getHeight()]; // x~column, y~row

		if (imagePlus.getType()== ImagePlus.GRAY8) {
			byte[] read = (byte[]) imagePlus.getProcessor().getPixels();
			for (int idx=0; idx<read.length; idx++) {
				inimg_xy[idx%imagePlus.getWidth()][idx/imagePlus.getWidth()] = (float) (read[idx] & 0xff);
			}

		}
		else if (imagePlus.getType()==ImagePlus.GRAY32) {
			float[] read = (float[]) imagePlus.getProcessor().getPixels();
			for (int idx=0; idx<read.length; idx++) {
				inimg_xy[idx%imagePlus.getWidth()][idx/imagePlus.getWidth()] = read[idx];
			}
		}
		else {
			IJ.log("image type not recognized");
			return DONE;
		}

		/******************************
        Generic Dialog
		 *****************************/
		nhoodRadius             = (float)   Prefs.get("advantra.critpoint.mask.nhoodRadius", 5);
		iDiff 					= (float)   Prefs.get("advantra.critpoint.mask.iDiff", 5);

		GenericDialog gd = new GenericDialog("MASKER2DDEMO");
		gd.addNumericField("radius ", 	nhoodRadius, 	0, 10, "spatial neighbourhood");
		gd.addNumericField("iDiff ", 	iDiff, 			0, 10, "intensity margin");

		gd.showDialog();
		if (gd.wasCanceled()) return DONE;

		nhoodRadius      = (float) gd.getNextNumber();
		Prefs.set("advantra.critpoint.mask.nhoodRadius",	nhoodRadius);

		iDiff       	= (float) gd.getNextNumber();
		Prefs.set("advantra.critpoint.mask.iDiff", 	    	iDiff);

		CPU_NR = Runtime.getRuntime().availableProcessors();

		return DOES_8G+DOES_32+NO_CHANGES;
    }

    public void run(ImageProcessor imageProcessor) {

        long t1, t2;

        /*
        main
         */
		t1 = System.currentTimeMillis();

        Masker2D.loadTemplate(inimg_xy, 0, nhoodRadius, iDiff);
        int totalLocs = inimg_xy.length * inimg_xy[0].length;

        Masker2D ms_jobs[] = new Masker2D[CPU_NR];
        for (int i = 0; i < ms_jobs.length; i++) {
            ms_jobs[i] = new Masker2D(i*totalLocs/CPU_NR,  (i+1)*totalLocs/CPU_NR);
            ms_jobs[i].start();
        }
        for (int i = 0; i < ms_jobs.length; i++) {
            try {
                ms_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

		Masker2D.formRemainingOutputs();

        t2 = System.currentTimeMillis();
        IJ.log("done. "+((t2-t1)/1000f)+"sec.");

        ImagePlus outmask = new ImagePlus("mask", Masker2D.getMask());
        outmask.show();

        ImagePlus outback = new ImagePlus("background", Masker2D.getBackground());
		outback.show();

		IJ.log("\ntotal "+Masker2D.i2xy.length+" locations extracted.\n"+"elapsed: "+((t2-t1)/1000f)+ " seconds.");

    }

//	public void mouseClicked(MouseEvent e) {
//		int offscreenX = incanvas.offScreenX(e.getX());
//		int offscreenY = incanvas.offScreenY(e.getY());
//		int radiusToCheck = (int) Math.ceil(nhoodRadius);
//		float[] circVals = Masker.extractCircleVals(offscreenX, offscreenY, radiusToCheck, Masker.circleElements(radiusToCheck));
//		Arrays.sort(circVals);
//		float[] x = new float[circVals.length];
//		for (int i=0; i<x.length; i++) x[i] = i+1;
//
//		Plot p = new Plot("", "", "", x, circVals);
//
//		float avg =  Masker.average(circVals);
//		float[] med = Masker.medianVec(circVals);
//		float[] std = Masker.stdVec(circVals, avg); for (int i=0; i<circVals.length; i++) std[i] = avg + 2*std[i];
//        float[] q3 = Masker.quartile3Vec(circVals);
//
//
//        float diff = std[0]-med[0];  // std , q3
//		float margin = (diff<=iDiff)?
//							   iDiff :
//							   0;//Masker.I_DIFF*(float)Math.exp(-0.5*(diff-Masker.I_DIFF)) ;
//		float[] mg = new float[circVals.length];
//		for (int i=0; i<circVals.length; i++) mg[i] = med[0] + margin;
//
//		p.setLimits(1, x.length, circVals[0], Math.max(circVals[circVals.length-1], Math.max(std[0], mg[0]))+5);
//
//        /*
//        BGRD
//         */
//
//		p.draw(); p.setLineWidth(2); p.setColor(Color.BLUE);
//		p.addPoints(x, med, Plot.LINE);
//
//        /*
//        LEVEL
//         */
//		p.draw(); p.setLineWidth(2); p.setColor(Color.GREEN);
//		p.addPoints(x, std, Plot.LINE);  // std, q3
//
//        /*
//        THRESHOLD
//         */
//		p.draw(); p.setLineWidth(2); p.setColor(Color.RED);
//		p.addPoints(x, mg, Plot.LINE);
//
//		IJ.log("median(R): "+med[0]+", upper limit(G): "+std[0]+", margin: "+mg[0]+", value:"+inimg.getProcessor().getf(offscreenX, offscreenY));
//
//		if (pw == null) pw = p.show();
//		pw.drawPlot(p);
//		pw.setTitle("x = " + offscreenX + ", y = " + offscreenY + ", r = " + nhoodRadius);
//
//		OvalRoi ring = new OvalRoi(offscreenX-nhoodRadius+.5, offscreenY-nhoodRadius+.5, 2*nhoodRadius, 2*nhoodRadius);
//		Overlay ov = new Overlay(); ov.add(ring);
//		incanvas.setOverlay(ov);
//
//	}
}
