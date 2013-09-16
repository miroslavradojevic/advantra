package detection;

import aux.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Overlay;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 9/10/13
 * Time: 7:57 PM
 */
public class JunctionDetection implements PlugInFilter {

	ImagePlus 	imp;
	ImageCanvas cnv;

	// parameters necessary for the detection
	double 	D;
	float 	iDiff;
	static int 	 	MIN_SIZE 	= 3;

	public int setup(String s, ImagePlus imagePlus) {

		if(imagePlus==null) return DONE;
		imp = Tools.convertToFloatImage(imagePlus);
		imp.setTitle("inimg");

		cnv = imagePlus.getCanvas();

		/***********************************************************/
		D       			=  			Prefs.get("advantra.critpoint.D", 3);
		iDiff				= (float) 	Prefs.get("advantra.critpoint.iDiff", 	5);

		GenericDialog gd = new GenericDialog("JUNCTION DET.");
		gd.addNumericField("neuronD ",  D, 		0, 5, " pix.");
		gd.addNumericField("iDiff ",  	iDiff, 	1, 5, " ");

		gd.showDialog();
		if (gd.wasCanceled()) return DONE;

		D       =  gd.getNextNumber();
		Prefs.set("advantra.critpoint.D",   D);

		iDiff               = (float) gd.getNextNumber();
		Prefs.set("advantra.critpoint.iDiff",   iDiff);
		/***********************************************************/

		return DOES_8G+DOES_32+NO_CHANGES;
	}

	public void run(ImageProcessor imageProcessor) {

		Detector det = new Detector();

		ArrayList<ArrayList<int[]>> detRegs = det.run(imp, D, iDiff);				// connected regions

		Vector<float[]> detLst = det.formDetectionList(detRegs, MIN_SIZE);			// list of discs

		int[] labels = det.clustering(detLst, Float.MAX_VALUE);  					// prune detection list, list of discs

		Vector<float[]> detLstPruned = Detector.groupClusters(labels, detLst);      // group clusters using labels

		Overlay detOv = det.formPointOverlay(detLstPruned, Color.RED);

		IJ.log("NEW, found " + detLstPruned.size());

		ImagePlus imNEW = new ImagePlus("NEW", imp.getProcessor());
		imNEW.setOverlay(detOv);
		imNEW.show();

		//***************************************************
		ArrayList<ArrayList<int[]>> detRegs_OLD = det.run_old(imp, D, iDiff);     	// connected regions

		Vector<float[]> detLst_OLD = det.formDetectionList(detRegs_OLD, MIN_SIZE);  // list of discs

		int[] labels_OLD = det.clustering(detLst_OLD, Float.MAX_VALUE);  							// prune detection list, list of discs

		Vector<float[]> detLstPruned_OLD = Detector.groupClusters(labels_OLD, detLst_OLD);// group clusters using labels

		Overlay detOv_OLD = det.formPointOverlay(detLstPruned_OLD, Color.YELLOW);

		IJ.log("OLD, found " + detLstPruned_OLD.size());

		ImagePlus imOLD = new ImagePlus("OLD", imp.getProcessor());
		imOLD.setOverlay(detOv_OLD);
		imOLD.show();


	}



}
