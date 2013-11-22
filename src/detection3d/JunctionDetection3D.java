package detection3d;

import ij.ImagePlus;
import ij.gui.ImageCanvas;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 11/22/13
 * Time: 6:03 AM
 */
public class JunctionDetection3D {

	ImagePlus imp;
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

	// in case it is necessary to downsample
	boolean downsample = false;
	float   neuriteDiameter = Float.NaN;



	// uses Detector3D class for the detection




}
