package detection3d;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Overlay;
import ij.gui.PointRoi;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 11/22/13
 * Time: 5:18 AM
 */
public class Detector3D {

    float[][][] img3d_zxy;  // inpput 3d image in convenient form (z,x,y)
    float zDist;

	Sphere3D sph3D;         // will be tool for thread classes - used to make geometric computations

    float 	D;
    float 	iDiff;
    float 	r;
    float 	scale = 1.5f;                   // heuristics: 1.5 is resonable

    public static float 	MIN_COS_ANG = .6f;
    public static float 	MIN_FUZZY_SCORE = .7f;
    public static float 	SCATTER_D2 	= 5;
    public static int		W_STD_RATIO_WRT_TO_D = 3;

    int CPU_NR;
//    float NBHOOD_SCALE = 1.0f;

    // static classes used for parallel processing
    Masker3D[]      masker_jobs     = null;
    MaskerOutput masker_output      = new MaskerOutput(); 	// foreground extraction is stored here
    Profiler3D[]    profiler_jobs   = null;
    PeakExtractor3D[] peak_extractor_jobs = null;

    public Detector3D() {
        System.out.println("loading default parameters...");
    }

	public Detector3D( // constructor initializes the parameter, run() calls different images with predefined parameters
            float   minCosAng,
            float   minFuzzyScore,
            float   scatterDistance2,
            int     wStdRatioToD
    ) {

		// store parameters
        MIN_COS_ANG = minCosAng;
        MIN_FUZZY_SCORE = minFuzzyScore;
        SCATTER_D2 = scatterDistance2;
        W_STD_RATIO_WRT_TO_D = wStdRatioToD;

    }

	public void run(ImagePlus inimg, float zDist1, float neuronDiameter, float iDiff1) { // ArrayList<ArrayList<int[]>>

        // load image into float[][][] structure
        if(inimg==null || inimg.getStackSize()<=1) return ;

        this.img3d_zxy = stackToZxyArray(inimg.getStack());
        this.zDist = zDist1 / 3;                                    // correct for the spreading through slices
        this.D = neuronDiameter;
        this.iDiff = iDiff1;
        this.r = scale * D;

		this.sph3D = new Sphere3D(D, scale);

        // some parameters
        this.CPU_NR = Runtime.getRuntime().availableProcessors();

        /*
            general pipeline is:

            1 find mask containing foreground voxels    (Masker3D)

            2 extract profiles in 3d for each location  (Profiler3D)

            3 find peaks of the profiles and correlate them to subpixel locations in 3d ()

            4 analyze recursive lists starting from each location
         */

        long t1, t2;




        //// 1. MASK
        System.out.println("separating foreground... ");
		t1 = System.currentTimeMillis();
        int margin = 0;
        Masker3D.loadTemplate(img3d_zxy, margin, sph3D.getOuterSamplingRadius()+2, zDist, iDiff);// (include shift 2 for the 4-neighbourhood check - to avoid out of array error)
        int totalMaskerJobs = img3d_zxy.length * img3d_zxy[0].length * img3d_zxy[0][0].length;
        masker_jobs = new Masker3D[CPU_NR];
        for (int i = 0; i < masker_jobs.length; i++) {
            masker_jobs[i] = new Masker3D(i*totalMaskerJobs/CPU_NR,  (i+1)*totalMaskerJobs/CPU_NR);
            masker_jobs[i].start();
        }
        for (int i = 0; i < masker_jobs.length; i++) {
            try {
                masker_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

		/*
		 	skipped locations are filled,
		*/
        Masker3D.fill();
        masker_output = Masker3D.getMaskerOutput();
        t2 = System.currentTimeMillis();
		System.out.println("done. " + ((t2 - t1) / 1000f) + " sec.");









		//// 2. PROFILE
        System.out.println("extracting profiles... ");
        t1 = System.currentTimeMillis();
        Profiler3D.loadTemplate(sph3D, masker_output.foregroundLocsZXY, img3d_zxy, zDist);
        int totalProfilerJobs = sph3D.getProfileLength();
        profiler_jobs = new Profiler3D[CPU_NR];
        for (int i = 0; i < profiler_jobs.length; i++) {
            profiler_jobs[i] = new Profiler3D(i*totalProfilerJobs/CPU_NR,  (i+1)*totalProfilerJobs/CPU_NR);
            profiler_jobs[i].start();
        }
        for (int i = 0; i < profiler_jobs.length; i++) {
            try {
                profiler_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
        t2 = System.currentTimeMillis();
        System.out.println("done. " + ((t2 - t1) / 1000f) + " sec.");








		//// 3. PROFILE PEAK(S)
        System.out.println("extracting peaks... ");
		System.out.println(""+img3d_zxy.length);
		System.out.println(""+img3d_zxy[0].length);
		System.out.println(""+img3d_zxy[0][0].length);

        t1 = System.currentTimeMillis();

        PeakExtractor3D.loadTemplate(sph3D, masker_output.foregroundLocsZXY, Profiler3D.prof3, img3d_zxy, zDist);
        int totalPeakExtractorJobs = masker_output.foregroundLocsZXY.length;
        peak_extractor_jobs = new PeakExtractor3D[CPU_NR];
        for (int i=0; i < peak_extractor_jobs.length; i++) {
			peak_extractor_jobs[i] = new PeakExtractor3D(i*totalPeakExtractorJobs/CPU_NR, (i+1)*totalPeakExtractorJobs/CPU_NR);
			peak_extractor_jobs[i].start();
        }
		for (int i=0; i < peak_extractor_jobs.length; i++) {
			try {
				peak_extractor_jobs[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
        t2 = System.currentTimeMillis();
        System.out.println("done. " + ((t2 - t1) / 1000f) + " sec.");










		int X = 261, Y = 137, Z = 27;
		int index = masker_output.locIndexZXY[Z][X][Y];
		System.out.println("index: "+index);
		System.out.println("check: Z: "+masker_output.foregroundLocsZXY[index][0]+", X: "+masker_output.foregroundLocsZXY[index][1]+", Y: "+masker_output.foregroundLocsZXY[index][2]);


		short[] curr_profile =  Profiler3D.prof3[index];
		ImagePlus img_profile = new ImagePlus("profile", sph3D.drawProfile(curr_profile));
		ArrayList<int[]> locs = sph3D.profilePeaksXY(curr_profile);
		Overlay ov_peaks = new Overlay();
		for (int ll = 0; ll<locs.size(); ll++) {
			ov_peaks.add(new PointRoi(locs.get(ll)[0]+.5, locs.get(ll)[1]+.5));
		}

		ArrayList<Integer> peakIdxs = sph3D.profilePeaks(curr_profile);
		for (int jj = 0; jj<peakIdxs.size(); jj++) { // masks.get(ii).length
			for (int ii=0; ii<sph3D.masks.get(peakIdxs.get(jj)).length; ii++) {

				int neighbourIdx = sph3D.masks.get(peakIdxs.get(jj))[ii];
				int neighbourXplot = sph3D.vizXY.get(neighbourIdx)[0];
				int neighbourYplot = sph3D.vizXY.get(neighbourIdx)[1];
//				layerMask.setf(neighbourXplot, neighbourYplot, 255f);
				PointRoi pt_mask = new PointRoi(neighbourXplot+.5, neighbourYplot+.5);
				pt_mask.setFillColor(Color.BLUE);
				ov_peaks.add(pt_mask);

			}
		}


		img_profile.setOverlay(ov_peaks);
		img_profile.show();


		Overlay ov = new Overlay();
		PointRoi pt = new PointRoi(X+0.5, Y+0.5);
		pt.setPosition(Z);
		pt.setStrokeColor(Color.RED);
		pt.setFillColor(Color.RED);
		ov.add(pt);

		for (int ii = 0; ii< 4; ii++) {
			System.out.println("peak "+ii+" -> "+Arrays.toString(PeakExtractor3D.peaks3[index][ii]));
			pt = new PointRoi(PeakExtractor3D.peaks3[index][ii][0]+.5, PeakExtractor3D.peaks3[index][ii][1]+.5);
			pt.setPosition(PeakExtractor3D.peaks3[index][ii][2]);
			ov.add(pt);
		}

		inimg.setOverlay(ov);
		inimg.show();

		new ImagePlus("masks", sph3D.visualizeMasks()).show();





        // 4. ASSOCIATE PEAKS



	}

    /*
        run threaded masker within the method
     */
//    private void runMasker(float[][][] img3d_zxy, float zDist, float radiusMask, float iDiff, MaskerOutput mo)
//    {
//		/*
//			following bit will execute run() in parallel
//		 */
////		int Zdim = img3d_zxy.length;
////		int Xdim = img3d_zxy[0].length;
////		int Ydim = img3d_zxy[0][0].length;
////        System.out.println("done with fill()");
//		/*
//		 	store locations, estimated background and binary mask them in output class (take from Masker3D into MaskerOutput instance
//		  */
////		mo.isForeground = new boolean[Zdim][Xdim][Ydim];
////		new ImagePlus("", Masker3D.getMask()).show();
//	}

	/*
		run threaded profile (peak really) extraction
	 */
//	private void runProfiler()
//	{
//
//		/*
//			run profiler in parallel - extract peaks
//		 */
//
//	}

	/*
		convert image to the array form (float[][][]) that will be used by all thread classes (Masker3D, Profiler3D, PeakExtractor, Analyzer3D) in run()
	 */
	private static float[][][] stackToZxyArray(ImageStack inis) {

		int W = inis.getWidth();
		int H = inis.getHeight();
		int L = inis.getSize();

		float[][][] img3d_zxy = new float[L][][];

		for (int l=0; l<L; l++) {

			img3d_zxy[l] = new float[W][H];
			float[] readSlice = (float[]) inis.getProcessor(l+1).convertToFloat().getPixels();

			for (int ww=0; ww<W; ww++) {
				for (int hh=0; hh<H; hh++) {
					img3d_zxy[l][ww][hh] = readSlice[hh*W+ww];
				}
			}

		}

		return img3d_zxy;
	}

    // TODO add method that converts opposite way float[][][] to ImageStack

    // TODO add method to export mask, background  as ImageStack, or ImagePlus

}