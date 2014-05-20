package detection2d;

import aux.ReadSWC;
import conn.Find_Connected_Regions;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.io.FileSaver;
import ij.process.ByteProcessor;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;

/**
 * Created by miroslav on 16-5-14.
 */
public class Detector2D {

	String 		image_dir;
	String		image_name;
	float[][] 	inimg_xy;               	// store input image as an array
	String 		image_gndtth_endpoints;		// ground truth swc file with critical points to evaluate, same name as input image
	String 		image_gndtth_bifurcations;




	// parameters
	float       	s;
	static float[]  D;
	int         	L;
	float       	ncc_high;
	float       	ncc_low;
	float 			likelihood_high;
	float 			likelihood_low;
	float			output_sigma;


	float       min_D;
	float       max_D;



	int			min_size_bif;                 	// smallest connected region to be valid critical point
	int			max_size_bif;                 	// smallest connected region to be valid critical point
	int			min_size_end;                 	// smallest connected region to be valid critical point
	int			max_size_end;                 	// smallest connected region to be valid critical point
	float 		output_membership_th;       	// based on k and output_sigma



	int         	M = 1;
	float       	minCos = -.5f;
	float       	sampling_crosswise = .3f;
	float			k = 0.5f;                   // hardcoded, defines the stricktness of the threshold for out membership function, k higher => means more strictness
	static  float   eval_scale  = 1.5f;
	String          eval_string = "";



	String		output_dir_name; 		    	// parameter coded output folder name
	String 		output_log_name;
	PrintWriter logWriter = null;


	int         CPU_NR;

	// output
	static float[][][] critpoint_det; 		// final off/on fuzzy memberships are stored and updated here  critpoint_det[0]-end, [1]-non, [2]-bif
	ArrayList<float[]> bif_list;			// list of bifurcation detections (x,y,size,avg_score)
	ArrayList<float[]> end_list;


    public Detector2D(
							 ImagePlus ip_load,
							 float _s,
							 float[] _D,
							 int _L,
							 float _ncc_high,
							 float _ncc_low,
							 float _likelihood_high,
							 float _likelihood_low,
							 float _output_sigma
    ){

		image_dir = ip_load.getOriginalFileInfo().directory; //  + File.separator  + image_name
		image_name = ip_load.getShortTitle();

		inimg_xy = new float[ip_load.getWidth()][ip_load.getHeight()]; // x~column, y~row
		if (ip_load.getType()== ImagePlus.GRAY8) {
			byte[] read = (byte[]) ip_load.getProcessor().getPixels();
			for (int idx=0; idx<read.length; idx++) {
				inimg_xy[idx%ip_load.getWidth()][idx/ip_load.getWidth()] = (float) (read[idx] & 0xff);
			}
		}
		else if (ip_load.getType()==ImagePlus.GRAY32) {
			float[] read = (float[]) ip_load.getProcessor().getPixels();
			for (int idx=0; idx<read.length; idx++) {
				inimg_xy[idx%ip_load.getWidth()][idx/ip_load.getWidth()] = read[idx];
			}
		}
		else {
			IJ.log("image type not recognized");
			return;
		}

		image_gndtth_endpoints = image_name + ".end";
		image_gndtth_bifurcations = image_name + ".bif";




		s = _s;
		D = _D; // just pointer assignement, no memory allocated
		L = _L;
		ncc_high = _ncc_high;
		ncc_low = _ncc_low;
		likelihood_high = _likelihood_high;
		likelihood_low = _likelihood_low;
		output_sigma = _output_sigma;





		min_D = Float.POSITIVE_INFINITY;
		max_D = Float.NEGATIVE_INFINITY;
		String		Dlist = ""; // for the out directory
		for (int i=0; i<_D.length; i++) {
			if (_D[i]<min_D) min_D = _D[i];
			if (_D[i]>max_D) max_D = _D[i];
			Dlist += IJ.d2s(_D[i],0) + ((i==_D.length-1)?"":":");
		}

		min_size_bif = Math.round(min_D);
		max_size_bif = Math.round(min_D*min_D);
		min_size_end = Math.round(min_D);
		max_size_end = Math.round(min_D * min_D);


		output_membership_th = (float) Math.exp(-(0.5f*0.5f)/(2*output_sigma*output_sigma)); // 0.5 is middle of the output range
		output_membership_th = (float) (1 - Math.pow(output_membership_th,k) * (1-output_membership_th)); // h1 = 1 - h^k * (1-h)



		output_dir_name = image_dir+String.format(       "det.s.Dlist.M.L.nncH.nccL.lhoodH.lhoodL.outSig._%.1f_%s_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f",
														 this.s,
														 Dlist,
														 M, L,
														 ncc_high,
														 ncc_low,
														 likelihood_high,
														 likelihood_low,
														 output_sigma);
		output_log_name = output_dir_name+ File.separator+"det.csv";

		File f = new File(output_dir_name);
		if (!f.exists()) {

			f.mkdirs();

			try {
				logWriter = new PrintWriter(output_log_name);
				logWriter.print("name,\tTP_BIF,\tFP_BIF,\tFN_BIF,\tP_BIF,\tR_BIF,\tTP_END,\tFP_END,\tFN_END,\tP_END,\tR_END\n");
				logWriter.close();
			} catch (FileNotFoundException ex) {}

		}
		// if it exists already in the folder, just prepare to append on the existing file
		try {
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(output_log_name, true)));
		} catch (IOException e) {}

		CPU_NR = Runtime.getRuntime().availableProcessors() + 1;

		// allocate outputs: off and on values for selected critpoint (static - to be updated throughout the process)
		critpoint_det = new float[3][ip_load.getWidth()][ip_load.getHeight()];

    }

    public void run(){

		long t1, t2;
		t1 = System.currentTimeMillis();

		for (int didx = 0; didx<D.length; didx++) {





			Sphere2D sph2d = new Sphere2D(D[didx], s); //sph2d.showSampling().show(); sph2d.showWeights().show();




			float new_masker_radius = 1.5f*sph2d.getOuterRadius();   // important that it is outer radius of the sphere
			float new_masker_percentile = 30;                   // used to have these two as argument but not necessary
			Masker2D.loadTemplate(inimg_xy, (int)Math.ceil(new_masker_radius), new_masker_radius, new_masker_percentile); //image, margin, check, percentile
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
			Masker2D.defineThreshold();
			Masker2D.formRemainingOutputs();
			//new ImagePlus("MASK,r="+IJ.d2s(new_masker_radius,2)+",per="+IJ.d2s(new_masker_percentile,0), Masker2D.getMask()).show();
			/********/









//			IJ.log("calculating profiles...");
			Profiler2D.loadTemplate(sph2d, Masker2D.i2xy, Masker2D.xy2i, inimg_xy);
			int totalProfileComponents = sph2d.getProfileLength();
			Profiler2D pf_jobs[] = new Profiler2D[CPU_NR];
			for (int i = 0; i < pf_jobs.length; i++) {
				pf_jobs[i] = new Profiler2D(i*totalProfileComponents/CPU_NR,  (i+1)*totalProfileComponents/CPU_NR);
				pf_jobs[i].start();
			}
			for (int i = 0; i < pf_jobs.length; i++) {
				try {
					pf_jobs[i].join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}





			/********/
//			IJ.log("extracting peaks...");
			PeakExtractor2D.loadTemplate(sph2d, Masker2D.i2xy, Profiler2D.prof2, inimg_xy, Masker2D.xy2i);
			int totalPeakExtrComponents = Profiler2D.prof2.length; // number of profiles == number of locations i2xy.length
			PeakExtractor2D pe_jobs[] = new PeakExtractor2D[CPU_NR];
			for (int i = 0; i < pe_jobs.length; i++) {
				pe_jobs[i] = new PeakExtractor2D(i*totalPeakExtrComponents/CPU_NR, (i+1)*totalPeakExtrComponents/CPU_NR);
				pe_jobs[i].start();
			}
			for (int i = 0; i < pe_jobs.length; i++) {
				try {
					pe_jobs[i].join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			//PeakExtractor2D.getCircStat().show();







			/********/
//			IJ.log("fitting model... extracting features... fuzzy detection...");
			FeatureExtractor2D.loadTemplate(
												   Masker2D.i2xy,
												   Masker2D.xy2i,
												   PeakExtractor2D.peaks_i,
												   PeakExtractor2D.peaks_w,
												   PeakExtractor2D.peaks_lhood,
												   inimg_xy,
												   M,
												   minCos,
												   D[didx],
												   "MEDIAN",
												   ncc_high,
												   ncc_low,
												   likelihood_high,
												   likelihood_low,
												   output_sigma,
												   L,
												   sampling_crosswise
			);
			int totalPeakAnalyzeComponents = Masker2D.i2xy.length; // number of locations
			FeatureExtractor2D pa_jobs[] = new FeatureExtractor2D[CPU_NR];
			for (int i = 0; i < pa_jobs.length; i++) {
				pa_jobs[i] = new FeatureExtractor2D(i*totalPeakAnalyzeComponents/CPU_NR, (i+1)*totalPeakAnalyzeComponents/CPU_NR);
				pa_jobs[i].start();
			}
			for (int i = 0; i < pa_jobs.length; i++) {
				try {
					pa_jobs[i].join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}





			/********/
//			IJ.log("append scores...");
			// take critpoint on/off scores from last class and append them to the critpoint_det[][][] - scale invariance
			for (int ll=0; ll<FeatureExtractor2D.critpoint_score2.length; ll++) {
				if (FeatureExtractor2D.critpoint_score2[ll]!=null) {

					float EE 		= FeatureExtractor2D.critpoint_score2[ll][0];
					float NN 	    = FeatureExtractor2D.critpoint_score2[ll][1];
					float BB 		= FeatureExtractor2D.critpoint_score2[ll][2];

					int x = Masker2D.i2xy[ll][0];  // masker will be redefined at each scale
					int y = Masker2D.i2xy[ll][1];

					float EEE 	= critpoint_det[0][x][y];
					float BBB 	= critpoint_det[2][x][y];

					if (EE > EEE || BB > BBB) { // it happened that one scale had higher score - take it
						critpoint_det[0][x][y] = EE;
						critpoint_det[1][x][y] = NN;
						critpoint_det[2][x][y] = BB;
					}
				}

			}



		}

		// extract the list of END and BIF locations

		bif_list = extractClusters(critpoint_det[2], min_size_bif, max_size_bif);
		end_list = extractClusters(critpoint_det[0], min_size_end, max_size_end);

		t2 = System.currentTimeMillis();
		IJ.log("done. " + ((t2 - t1) / 1000f) + "sec.");
		IJ.log(bif_list.size()+" bifurcations ");
		IJ.log(end_list.size()+" endpoints ");

	}

	public void doEvaluation(){

		eval_string = ""; // initialize before appending in following calls

		int tp, fp, fn;
		float prec, recl;
		ReadSWC reader;
		boolean[] annots;

		// evaluate bifurcations
		if (new File(image_dir+image_gndtth_bifurcations).exists()) {

//            IJ.log("evaluating BIFURCATION detection using " + image_dir+image_gndtth_bifurcations);
			tp = 0; fp = 0; fn = 0;
			reader = new ReadSWC(image_dir+image_gndtth_bifurcations);
			annots = new boolean[reader.nodes.size()];

			for (int a=0; a<bif_list.size(); a++) {

				boolean found = false;

				for (int b=0; b<reader.nodes.size(); b++) {
					float ax = bif_list.get(a)[0];
					float ay = bif_list.get(a)[1];
					float bx = reader.nodes.get(b)[reader.XCOORD];
					float by = reader.nodes.get(b)[reader.YCOORD];
					float dist = (float) Math.sqrt(Math.pow(bx-ax,2)+Math.pow(by-ay,2));
					if (dist<=eval_scale*max_D) {// compare the distance
						found = true;
						tp++;
						annots[b] = true;
					}
				}

				if (!found) fp++;

			}

			for (int a=0; a<annots.length; a++) if (!annots[a]) fn++;

			prec = (tp+fp>0)? (tp/(float)(tp+fp)) : 0 ;
			recl = (tp+fn>0)? (tp/(float)(tp+fn)) : 0 ;
			eval_string += "\""+image_name+"\",\t" + tp + ",\t" + fp + ",\t"+fn+ ",\t" + IJ.d2s(prec, 2) + ",\t" + IJ.d2s(recl,2) +",\t";

//            IJ.log("evaluating ENDPOINT detection using " + image_dir+image_gndtth_endpoints);
			tp = 0; fp = 0; fn = 0;
			reader = new ReadSWC(image_dir+image_gndtth_endpoints);
			annots = new boolean[reader.nodes.size()];

			for (int a=0; a<end_list.size(); a++) {

				boolean found = false;

				for (int b=0; b<reader.nodes.size(); b++) {
					float ax = end_list.get(a)[0];
					float ay = end_list.get(a)[1];
					float bx = reader.nodes.get(b)[reader.XCOORD];
					float by = reader.nodes.get(b)[reader.YCOORD];
					float dist = (float) Math.sqrt(Math.pow(bx-ax,2)+Math.pow(by-ay,2));
					if (dist<=eval_scale*max_D) {// compare the distance
						found = true;
						tp++;
						annots[b] = true;
					}
				}

				if (!found) fp++;

			}

			for (int a=0; a<annots.length; a++) if (!annots[a]) fn++;

			prec = (tp+fp>0)? (tp/(float)(tp+fp)) : 0 ;
			recl = (tp+fn>0)? (tp/(float)(tp+fn)) : 0 ;
			eval_string += tp + ",\t" + fp + ",\t"+fn+ ",\t" + IJ.d2s(prec,2) + ",\t" + IJ.d2s(recl,2);

		}
		else
			eval_string += "\"" + image_name + "\",\t" + Float.NaN + ",\t" + Float.NaN + ",\t" + Float.NaN + ",\t" + Float.NaN + ",\t" + Float.NaN + ",\t" +
								   						 Float.NaN + ",\t" + Float.NaN + ",\t" + Float.NaN + ",\t" + Float.NaN + ",\t" + Float.NaN;

		logWriter.println(eval_string);
		IJ.log(""+eval_string);

		logWriter.close();

	}

	public Overlay getBifurcationOverlay()
	{

		// create red Overlay with BIFs
		Overlay ov = new Overlay();
		for (int i=0; i<bif_list.size(); i++) {
			// x,y,size,avg_score in bif_list
			float cx = bif_list.get(i)[0];
			float cy = bif_list.get(i)[1];
			float sc = bif_list.get(i)[3];

			float R = 0.5f * max_D;
			OvalRoi ovroi = new OvalRoi(cx-R+.5, cy-R+.5, 2*R, 2*R);
			ovroi.setStrokeWidth(3);
			ovroi.setFillColor(new Color(1, 0, 0, sc));
			ov.add(ovroi);
		}

		return ov;

	}

	public Overlay getEndpointOverlay()
	{

        // create yellow Overlay with ENDs
		Overlay ov = new Overlay();//formPointOverlay(conn_reg.getConnectedRegions(), min_size_end, max_size_end, 1, 1, 0); // add yellow intensity based on average score
		for (int i=0; i<end_list.size(); i++) {
			// x,y,size,avg_score in bif_list
			float cx = end_list.get(i)[0];
			float cy = end_list.get(i)[1];
			float sc = end_list.get(i)[3];

			float R = 0.5f * max_D;
			OvalRoi ovroi = new OvalRoi(cx-R+.5, cy-R+.5, 2*R, 2*R);
			ovroi.setStrokeWidth(3);
			ovroi.setFillColor(new Color(1, 1, 0, sc));
			ov.add(ovroi);
		}

		return ov;

	}

	private ArrayList<float[]> extractClusters(float[][] _critpoint_det, int minSize, int maxSize)
	{// x,y,size,avg_score

		// create detection image
		int w = _critpoint_det.length;
		int h = _critpoint_det[0].length;

		byte[] t = new byte[w*h];
		for (int x=0; x<w; x++) {
			for (int y=0; y<h; y++) {
				float curr_on = _critpoint_det[x][y];
				if (curr_on>=output_membership_th) {
					t[y*w+x] = (byte) 255;
				}
			}
		}

		// take detections (binary image), find connected regions
		ByteProcessor bp = new ByteProcessor(w, h, t);
		Find_Connected_Regions conn_reg = new Find_Connected_Regions(new ImagePlus("", bp), true);
		conn_reg.run("");
		//conn_reg.showLabels().show();
		//FileSaver fs = new FileSaver(conn_reg.showLabels());
		//fs.saveAsTiff(image_dir+"conn_ends.tif");
		ArrayList<ArrayList<int[]>> regs = conn_reg.getConnectedRegions();
		ArrayList<float[]> cls = new ArrayList<float[]>(); // x,y,size,avg_score

		// regs are converted into cls
		for (int i=0; i<regs.size(); i++) {
			if (regs.get(i).size()>=minSize && regs.get(i).size()<=maxSize) {

				float Cx=0, Cy=0;
				float C=0;

				for (int aa=0; aa<regs.get(i).size(); aa++) {

					int xcoord = regs.get(i).get(aa)[1];
					int ycoord = regs.get(i).get(aa)[0];

					Cx += regs.get(i).get(aa)[1];
					Cy += regs.get(i).get(aa)[0];
					C += _critpoint_det[xcoord][ycoord];

				}

				Cx /= regs.get(i).size();
				Cy /= regs.get(i).size();
				C /= regs.get(i).size();

				cls.add(new float[]{Cx, Cy, regs.get(i).size(), C});

			}
		}

		return cls;

	}

}