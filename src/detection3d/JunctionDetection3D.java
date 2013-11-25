package detection3d;

import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 11/22/13
 * Time: 6:03 AM
 */
public class JunctionDetection3D implements PlugInFilter {

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

    public int setup(String s, ImagePlus imagePlus) {
        return 0;
    }

    public void run(ImageProcessor imageProcessor) {
        //
    }

    public static void main(String[] args){

        ImagePlus inimg  = new ImagePlus("/home/miroslav/mosaic.tif");
        float zDist     = 3.4f;
        float D         = 5;
        float iDiff     = 10;

        // detector parameters

        // uses Detector3D
        Detector3D det3D = new Detector3D();
        det3D.run(inimg, zDist, D, iDiff);


    }

}
