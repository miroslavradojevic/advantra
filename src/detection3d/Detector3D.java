package detection3d;

import conn.Find_Connected_Regions;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.process.ShortProcessor;

import java.awt.*;
import java.io.*;
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
	Fuzzy3D fz3d; 			// will be tool for detection

    float 	D;
    float 	iDiff;
    float 	r;
    float 	scale = 1.5f;                   // heuristics: 1.5 is resonable

    public static float 	MIN_COS_ANG = .6f;
    public static float 	MIN_FUZZY_SCORE = .5f;
    public static float 	SCATTER_D2 	= 5;
    public static int		W_STD_RATIO_WRT_TO_D = 3;

	public static int 		MIN_SIZE = 1;

    int CPU_NR;
//    float NBHOOD_SCALE = 1.0f;

    // static classes used for parallel processing
    Masker3D[]      masker_jobs     = null;
    MaskerOutput masker_output      = new MaskerOutput(); 	// foreground extraction is stored here
    Profiler3D[]    profiler_jobs   = null;
    PeakExtractor3D[] peak_extractor_jobs = null;
    PeakAnalyzer3D[] peak_analyzer_jobs = null;
    ScoreCalculator3D[] score_calc_jobs = null;

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
	                //ArrayList<float[]>
	public float[][] run(ImagePlus inimg, float zDist1, float neuronDiameter, float iDiff1) { // ArrayList<ArrayList<int[]>>    // ArrayList<ArrayList<int[]>>

        // load image into float[][][] structure

		if(inimg==null || inimg.getStackSize()<=1) return null;

        this.img3d_zxy = stackToZxyArray(inimg.getStack());
        this.zDist = zDist1 / 5;                                    // correct for the spreading through slices
        this.D = neuronDiameter;
        this.iDiff = iDiff1;
        this.r = scale * D;

		this.sph3D = new Sphere3D(D, scale);
		this.fz3d = new Fuzzy3D(iDiff);

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

        /*
        ********************************************************
        * EXTRACT FOREGROUND
        ********************************************************
         */
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
        Masker3D.fill();  // skipped locations are filled
        masker_output = Masker3D.getMaskerOutput();
        t2 = System.currentTimeMillis();
		System.out.println("done. " + ((t2 - t1) / 1000f) + " sec.");
        //new ImagePlus("FG", Masker3D.getMask()).show();

		//exportSwcForegroundLocations("/home/miroslav/fgr.swc");
		//if (true) {System.out.println("finishinf"); return null;}

		/*
        ********************************************************
        * PROFILE EXTRACT
        ********************************************************
         */
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
        /*
        ********************************************************
        * PROFILE PEAK(S)
        ********************************************************
         */
        System.out.println("extracting peaks... ");
        t1 = System.currentTimeMillis();
        PeakExtractor3D.loadTemplate(sph3D, masker_output.foregroundLocsZXY, Profiler3D.prof3, img3d_zxy, zDist, masker_output.locIndexZXY);
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
        /*
        ********************************************************
        * PEAK(S) ANALYZE
        ********************************************************
         */
        System.out.println("associating peaks... ");
        t1 = System.currentTimeMillis();
        PeakAnalyzer3D.loadTemplate(sph3D, masker_output.foregroundLocsZXY, PeakExtractor3D.peaks3, masker_output.locIndexZXY, img3d_zxy, MIN_COS_ANG);
        int totalPeakAnalyzerJobs = masker_output.foregroundLocsZXY.length;
        peak_analyzer_jobs = new PeakAnalyzer3D[CPU_NR];
        for (int i=0; i < peak_analyzer_jobs.length; i++) {
            peak_analyzer_jobs[i] = new PeakAnalyzer3D(i*totalPeakAnalyzerJobs/CPU_NR, (i+1)*totalPeakAnalyzerJobs/CPU_NR);
            peak_analyzer_jobs[i].start();
        }
        for (int i=0; i < peak_analyzer_jobs.length; i++) {
            try {
                peak_analyzer_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
        t2 = System.currentTimeMillis();
        System.out.println("done. " + ((t2 - t1) / 1000f) + " sec.");
        /*
        ********************************************************
        * SCORE CALCULATION
        ********************************************************
         */
        System.out.println("calculating scores... ");
        t1 = System.currentTimeMillis();
        ScoreCalculator3D.loadTemplate(PeakAnalyzer3D.delin3, img3d_zxy, masker_output.locIndexZXY, masker_output.foregroundLocsZXY, Masker3D.back3, fz3d);
        int totalScoreCalcJobs = PeakAnalyzer3D.delin3.length;
        score_calc_jobs = new ScoreCalculator3D[CPU_NR];
        for (int i=0; i < score_calc_jobs.length; i++) {
            score_calc_jobs[i] = new ScoreCalculator3D(i*totalScoreCalcJobs/CPU_NR, (i+1)*totalScoreCalcJobs/CPU_NR);
            score_calc_jobs[i].start();
        }
        for (int i=0; i < score_calc_jobs.length; i++) {
            try {
                score_calc_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
		t2 = System.currentTimeMillis();
        System.out.println("done. " + ((t2 - t1) / 1000f) + " sec.");
		/*
        ********************************************************
        * CONNECTED COMPONENTS
        ********************************************************
         */
		System.out.println("connected components... ");
		ImagePlus thded_scores = ScoreCalculator3D.getThresholdedScores(MIN_FUZZY_SCORE);
		thded_scores.show();
		Find_Connected_Regions conn_reg = new Find_Connected_Regions(thded_scores, true);  // true means save locations
		conn_reg.run("");
		ArrayList<ArrayList<int[]>> dets = conn_reg.getConnectedRegions3D_XYZ();
		System.out.println(" done. "+dets.size()+ " elements:");

		for (int k=0; k< dets.size(); k++) {
			System.out.println(dets.get(k).size()+" el.");
			for (int kk = 0; kk<dets.get(k).size();kk++) {
				System.out.println(Arrays.toString(dets.get(k).get(kk)));
			}
		}

		float[][] out = formDetectionListXYZR(dets, MIN_SIZE);
		return out;

//		System.out.println(a.size()+" regions");
		//fz3d.decisionCube().show();
//		System.out.println("DONE ALL");
        //new ImagePlus("", sph3D.visualizeMasks()).show();

	}


	// debugging

//	public Overlay exportOverlayForegroundLocations() {
//
//		Overlay ov = new Overlay();
//
//		for (int i=0; i<masker_output.foregroundLocsZXY.length; i++) {
//
//			int x = masker_output.foregroundLocsZXY[i][1];
//			int y = masker_output.foregroundLocsZXY[i][2];
//			int z = masker_output.foregroundLocsZXY[i][0];
//
//			PointRoi p = new PointRoi(x+0.5, y+0.5);
//			p.setPosition(z+1);
//			ov.add(p);
//
//		}
//
//		return ov;
//	}

	public Overlay exportOverlayPeakLocations(){

		Overlay ov = new Overlay();

		for (int i=0; i<PeakAnalyzer3D.delin3.length; i++) {

			int cnt = 0;

			for (int j=0; j<PeakAnalyzer3D.delin3[i].length; j++){

				boolean complete = true;
				for (int k = 0; k<PeakAnalyzer3D.delin3[i][j].length; k++ ){
					if (PeakAnalyzer3D.delin3[i][j][k]==-1){
						complete = false;
					}
				}

				if (complete) cnt++;

//				if(PeakExtractor3D.peaks3[i][j][0]!=-1){
//					cnt++;
//				}


			}

			if (cnt >= 3) {

				int x = PeakExtractor3D.listLocs3D[i][1]; //peaks3[i][j][0];
				int y = PeakExtractor3D.listLocs3D[i][2]; //peaks3[i][j][1];
				int z = PeakExtractor3D.listLocs3D[i][0]; //peaks3[i][j][2];

				PointRoi p = new PointRoi(x+0.5, y+0.5);
				p.setPosition(z+1);
				ov.add(p);
			}

		}

		return ov;
	}




	public Overlay getLocalSkeleton(int atX, int atY, int atZ) {  // make it swc export

		Overlay ov = new Overlay();
        float plot_radius = 0.5f;

        int locationIndex = masker_output.locIndexZXY[atZ][atX][atY];

		if (locationIndex!=-1) {

            // central point exists for sure
            int cx = atX;
            int cy = atY;
            int cz = atZ;
            OvalRoi ovroi = new OvalRoi(cx-plot_radius+.5, cy-plot_radius+.5, 2*plot_radius, 2*plot_radius);
            ovroi.setPosition(cz);
            ovroi.setFillColor(Color.RED);
            ov.add(ovroi);

			// check all that were different from -1
			int[][] skel = PeakAnalyzer3D.delin3[locationIndex];

            //IJ.log("\n");
            //for (int k=0; k<PeakAnalyzer3D.delin3[locationIndex].length; k++) {
            //    IJ.log( (k+1)+"thread -> "+Arrays.toString( PeakAnalyzer3D.delin3[locationIndex][k] ));
            //}
            //IJ.log("\n");

            for (int thr_id = 0; thr_id<skel.length; thr_id++) {

                for (int pt_id = 0; pt_id<skel[thr_id].length; pt_id++) {

                    if (skel[thr_id][pt_id]==-1) {
                        //break;
                    }
                    else {

                        cx = masker_output.foregroundLocsZXY[skel[thr_id][pt_id]][1];
                        cy = masker_output.foregroundLocsZXY[skel[thr_id][pt_id]][2];
                        cz = masker_output.foregroundLocsZXY[skel[thr_id][pt_id]][0];
                        ovroi = new OvalRoi(cx-plot_radius+.5, cy-plot_radius+.5, 2*plot_radius, 2*plot_radius);
                        ovroi.setPosition(cz);
                        ov.add(ovroi);

                    }

                }

            }

		}

        return ov;

	}

	public ShortProcessor getLocalProfile(int atX, int atY, int atZ) {  // a wrapper for Sphere3D method drawProfile()

		if (masker_output.locIndexZXY[atZ][atX][atY]!=-1) {

			return sph3D.drawProfile(Profiler3D.prof3[masker_output.locIndexZXY[atZ][atX][atY]]);

		}
		else {
			return null;// return the empty profile, we're in background
		}

	}

    public Overlay getLocalProfilePeaks(int atX, int atY, int atZ) {  // wrapper for Sphere3D method

        // will give all peaks

        Overlay ov = new Overlay();
        float radius_plot = .4f;

        if (masker_output.locIndexZXY[atZ][atX][atY]!=-1) {

            ArrayList<int[]> pks = sph3D.profilePeaksXY(Profiler3D.prof3[masker_output.locIndexZXY[atZ][atX][atY]]);

            for (int a=0; a<pks.size(); a++) {
                OvalRoi ovroi = new OvalRoi(pks.get(a)[0]-radius_plot+.5, pks.get(a)[1]-radius_plot+.5, 2*radius_plot, 2*radius_plot);
                ovroi.setFillColor(Color.RED);
                ov.add(ovroi);
            }

        }

        return ov;

    }

    public Overlay getSelectedLocalProfilePeaks(int atX, int atY, int atZ) {

        // will give selected peaks

        Overlay ov = new Overlay();
        float radius_plot = .75f;

        int locationID = masker_output.locIndexZXY[atZ][atX][atY];

        if (locationID!=-1) {

            for (int a = 0; a<PeakExtractor3D.peaks2[locationID].length; a++) {

                int cx = PeakExtractor3D.peaks2[locationID][a][0];
                int cy = PeakExtractor3D.peaks2[locationID][a][1];

                OvalRoi ovroi = new OvalRoi(cx-radius_plot+.5, cy-radius_plot+.5, 2*radius_plot, 2*radius_plot);
                ov.add(ovroi);

            }

        }

        return ov;

    }

    public void debug(int atX, int atY, int atZ){

        int locationID = masker_output.locIndexZXY[atZ][atX][atY];

        if (locationID!=-1) {

            PeakExtractor3D.summary(atX, atY, atZ);

            PeakAnalyzer3D.summary(atX, atY, atZ);

			// TODO plot fls output

        }

    }

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

	private static float[][] formDetectionListXYZR(ArrayList<ArrayList<int[]>> regs, int minSize){  // regs is considered as list of [xyz] // , float[][] out

		int cnt = 0;
		for (int i=0; i<regs.size(); i++) {
			if (regs.get(i).size()>=minSize)
				cnt++;
		}

		if (cnt<=0) return null;

		float[][] out = new float[cnt][4];


		cnt =0;
		for (int i=0; i<regs.size(); i++) {
			if (regs.get(i).size()>=minSize) {

				float volume = (float)regs.get(i).size();
				float Cx=0, Cy=0, Cz = 0;
				float R = (float) Math.pow((3*volume)/(4*Math.PI), 1f/3);
				R = (R<1)? 1 : R ;

				for (int aa=0; aa<regs.get(i).size(); aa++) {
					Cx += regs.get(i).get(aa)[0];   // they are stored as x,y.z in regs
					Cy += regs.get(i).get(aa)[1];
					Cz += regs.get(i).get(aa)[2];
				}

				Cx /= regs.get(i).size();
				Cy /= regs.get(i).size();
				Cz /= regs.get(i).size();

				//out.add(new float[]{Cx, Cy, Cz, R});
				out[cnt][0] = Cx;
				out[cnt][1] = Cy;
				out[cnt][2] = Cz;
				out[cnt][3] = R;

				cnt++;

			}
		}

		return out;

	}

	public void exportSwcForegroundLocations(String file_path){

		PrintWriter logWriter = null;

		try {
			logWriter = new PrintWriter(file_path);
			logWriter.print("");
			logWriter.close();
		}
		catch (FileNotFoundException ex) {}

		try {
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(file_path, true)));
			logWriter.println("# SWC WITH FG");
		} catch (IOException e) {}


		for (int i=0; i<masker_output.foregroundLocsZXY.length; i++) {

			int x = masker_output.foregroundLocsZXY[i][1];
			int y = masker_output.foregroundLocsZXY[i][2];
			int z = masker_output.foregroundLocsZXY[i][0];

			logWriter.println((i+1)+" "+0+" "+x+" "+y+" "+z+" "+0.5+" "+(-1));

		}

		logWriter.close();

	}

    // TODO add method that converts opposite way float[][][] to ImageStack
	// TODO add the tools that would prune the list and export swc with detections

}