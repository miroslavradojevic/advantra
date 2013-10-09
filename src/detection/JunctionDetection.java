package detection;

import aux.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Overlay;
import ij.measure.ResultsTable;
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

	float 	minCosAngle;
    float 	minFuzzyScore;
    int 	MIN_SIZE;
    float 	scatterDistSquared = 5f;

	static int 			wStdRatioToD		= 6;
    static float 	 	LOCATION_TOLERANCE_SCALE 	= 1.8f;

    private static float 	Deg2Rad = (float) (Math.PI/180f);
    private static float 	Rad2Deg = (float) (180f/Math.PI);

	public int setup(String s, ImagePlus imagePlus) {

		if(imagePlus==null) return DONE;
		imp = Tools.convertToFloatImage(imagePlus);
		imp.setTitle("inimg");

		cnv = imagePlus.getCanvas();

        String  MENU_minCosAngle = "";
        String[] minCosAngleOptions = new String[6];
        minCosAngleOptions[0] = "60";
        minCosAngleOptions[1] = "50";
        minCosAngleOptions[2] = "40";
        minCosAngleOptions[3] = "30";
        minCosAngleOptions[4] = "20";
        minCosAngleOptions[5] = "10";

        String  MENU_minFuzzyScore = "";
        String[] minFuzzyScoreOptions = new String[3];
        minFuzzyScoreOptions[0] = "0.6";
        minFuzzyScoreOptions[1] = "0.7";
        minFuzzyScoreOptions[2] = "0.8";

        String   MENU_min_size = "";
        String[] min_sizeOptions = new String[5];
        min_sizeOptions[0] = "1";
        min_sizeOptions[1] = "2";
        min_sizeOptions[2] = "3";
        min_sizeOptions[3] = "4";
        min_sizeOptions[4] = "5";

        /***********************************************************/
		D       			= (float)	Prefs.get("advantra.critpoint.D", 3);
		iDiff				= (float) 	Prefs.get("advantra.critpoint.iDiff", 	5);
        MENU_minCosAngle    =           Prefs.get("advantra.critpoint.minCosAngle", minCosAngleOptions[0]);
        MENU_minFuzzyScore  =           Prefs.get("advantra.critpoint.minFuzzyScore", minFuzzyScoreOptions[0]);
        MENU_min_size       =           Prefs.get("advantra.critpoint.MIN_SIZE", min_sizeOptions[0]);

		GenericDialog gd = new GenericDialog("JUNCTION DET.");
		gd.addNumericField("neuronD ",  D, 		0, 5, " pix.");
		gd.addNumericField("iDiff ",  	iDiff, 	1, 5, " ");
		gd.addChoice("min cos ",    minCosAngleOptions,     MENU_minCosAngle);
		gd.addChoice("min fuzzy ",  minFuzzyScoreOptions,   MENU_minFuzzyScore);
		gd.addChoice("min size ",   min_sizeOptions,        MENU_min_size);

		gd.showDialog();
		if (gd.wasCanceled()) return DONE;

		D       = (float) gd.getNextNumber();
		Prefs.set("advantra.critpoint.D",   D);

		iDiff               = (float) gd.getNextNumber();
		Prefs.set("advantra.critpoint.iDiff",   iDiff);

        MENU_minCosAngle = gd.getNextChoice();
        minCosAngle = (float) Math.cos(Float.valueOf(MENU_minCosAngle) * Deg2Rad);
        Prefs.set("advantra.critpoint.minCosAngle", MENU_minCosAngle);

        MENU_minFuzzyScore = gd.getNextChoice();
        minFuzzyScore = Float.valueOf(MENU_minFuzzyScore);
        Prefs.set("advantra.critpoint.minFuzzyScore", MENU_minFuzzyScore);

        MENU_min_size = gd.getNextChoice();
        MIN_SIZE = Integer.valueOf(MENU_min_size);
        Prefs.set("advantra.critpoint.MIN_SIZE", MENU_min_size);

        IJ.log("loading parameters...");
        IJ.log("neuron diameter \t\t"+D);
        IJ.log("intensity difference \t\t"+iDiff);
        IJ.log("max. divergence angle [deg] \t\t"+MENU_minCosAngle+", "+"cos(angle) \t\t"+minCosAngle);
        IJ.log("min fuzzy score \t\t"+minFuzzyScore);
        IJ.log("min size \t\t"+MIN_SIZE);
        IJ.log("-------------------------------");
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

        ResultsTable resTab = det.formResultsTable(detRegs, det.fuzzyScores, MIN_SIZE);

        resTab.show("BIFURCATIONS");

        det.formChimeraScript(detRegs, det.fuzzyScores, MIN_SIZE);

		//System.out.println("### "+detLst.size() + "detections, " + detLstPruned.size() + "pruned detections");
		//IJ.log("NEW, found " + detLstPruned.size());

		ImagePlus imNEW = new ImagePlus("DET", imp.getProcessor());
		imNEW.setOverlay(detOv);
		imNEW.show();

		//***************************************************
        if (false) {
            ArrayList<ArrayList<int[]>> detRegs_OLD = det.run_old(imp, D, iDiff);     	// connected regions
            Vector<float[]> detLst_OLD = det.formDetectionList(detRegs_OLD, MIN_SIZE);  // list of discs
            int[] labels_OLD = det.clustering(detLst_OLD, LOCATION_TOLERANCE_SCALE*D);  // prune detection list, list of discs
            Vector<float[]> detLstPruned_OLD = Detector.groupClusters(labels_OLD, detLst_OLD);// group clusters using labels
            Overlay detOv_OLD = det.formPointOverlay(detLstPruned_OLD, Color.YELLOW);
            IJ.log("OLD, found " + detLstPruned_OLD.size());
            ImagePlus imOLD = new ImagePlus("OLD", imp.getProcessor());
            imOLD.setOverlay(detOv_OLD);
            imOLD.show();
        }
        //***************************************************

	}

}
