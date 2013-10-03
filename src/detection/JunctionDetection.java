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
	float 	D;
	float 	iDiff;

	// TODO make these params
	static float 		minCosAngle 		= .8f;
	static float 		minFuzzyScore 		= .6f;
	static float 		scatterDistSquared 	=  5f;

	static int 			wStdRatioToD		= 3;
	static int 	 		MIN_SIZE 	= 1;

	static float 	 	LOCATION_TOLERANCE_SCALE 	= 1.5f;

	public int setup(String s, ImagePlus imagePlus) {

		if(imagePlus==null) return DONE;
		imp = Tools.convertToFloatImage(imagePlus);
		imp.setTitle("inimg");

		cnv = imagePlus.getCanvas();

		/***********************************************************/
		D       			= (float)	Prefs.get("advantra.critpoint.D", 3);
		iDiff				= (float) 	Prefs.get("advantra.critpoint.iDiff", 	5);

		GenericDialog gd = new GenericDialog("JUNCTION DET.");
		gd.addNumericField("neuronD ",  D, 		0, 5, " pix.");
		gd.addNumericField("iDiff ",  	iDiff, 	1, 5, " ");

		gd.showDialog();
		if (gd.wasCanceled()) return DONE;

		D       = (float) gd.getNextNumber();
		Prefs.set("advantra.critpoint.D",   D);

		iDiff               = (float) gd.getNextNumber();
		Prefs.set("advantra.critpoint.iDiff",   iDiff);
		/***********************************************************/

		return DOES_8G+DOES_32+NO_CHANGES;
	}

	public void run(ImageProcessor imageProcessor) {

		//Detector det = new Detector();
		Detector det = new Detector(minCosAngle, minFuzzyScore, scatterDistSquared, wStdRatioToD);

		ArrayList<ArrayList<int[]>> detRegs = det.run(imp, D, iDiff);				// connected regions

		Vector<float[]> detLst = det.formDetectionList(detRegs, MIN_SIZE);			// list of discs

		int[] labels = det.clustering(detLst, LOCATION_TOLERANCE_SCALE*D);  		// prune detection list, list of discs

		Vector<float[]> detLstPruned = Detector.groupClusters(labels, detLst);      // group clusters using labels

		Overlay detOv = det.formPointOverlay(detLstPruned, Color.RED);

		//System.out.println("### "+detLst.size() + "detections, " + detLstPruned.size() + "pruned detections");
		IJ.log("NEW, found " + detLstPruned.size());

		ImagePlus imNEW = new ImagePlus("NEW", imp.getProcessor());
		imNEW.setOverlay(detOv);
		imNEW.show();

		//***************************************************
		ArrayList<ArrayList<int[]>> detRegs_OLD = det.run_old(imp, D, iDiff);     	// connected regions

		Vector<float[]> detLst_OLD = det.formDetectionList(detRegs_OLD, MIN_SIZE);  // list of discs

		int[] labels_OLD = det.clustering(detLst_OLD, LOCATION_TOLERANCE_SCALE*D);  							// prune detection list, list of discs

		Vector<float[]> detLstPruned_OLD = Detector.groupClusters(labels_OLD, detLst_OLD);// group clusters using labels

		Overlay detOv_OLD = det.formPointOverlay(detLstPruned_OLD, Color.YELLOW);

		IJ.log("OLD, found " + detLstPruned_OLD.size());

		ImagePlus imOLD = new ImagePlus("OLD", imp.getProcessor());
		imOLD.setOverlay(detOv_OLD);
		imOLD.show();


	}



}
