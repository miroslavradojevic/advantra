package generate;

import aux.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.io.FileSaver;
import ij.process.FloatProcessor;
import imagescience.random.PoissonGenerator;

import java.awt.*;
import java.io.*;
import java.util.Random;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 9/6/13
 * Time: 1:13 PM
 */
public class Generator {

	/*
		generates:
		- image stack with grid of bifurcation points
		 bifurcations can have different SNR, branch diameters (D), variation of branch diameteres (dD)
		(option: points can be also randomly distributed - but it's not useful for what I need them for)
		different angular configurations and branch diameters, gaussian shaped cross-section intensity profile
		(option: convolved with PSF - NOT NECESSARY really) and
		added Poisson noise
		- stores the ground truth locations in csv file
		- image with collection of bifurcation points (SNR, branch diameters as parameteres)
		- stores the ground truth locations in csv file
	 */

	// pixel width  = 0.3125849 um, or 84.66 um, two examples I got
	// NA = 1.4
	// lambda = 600 nm (525 actually)

	// to emulate PSF effect: IJ.run(imp, "Gaussian Blur...", "sigma=1.48");
	// Gaussian blur standard deviation:

	// 1. generate binary structure (b, h)
	// 2. PSF (Gaussian blur with appropriate sigma ~ 0.21 * lambda / NA) => low res sigmaXY = 9.3e-4, high res 0.252
	// 3. generate Poisson Noise wrt to the image values

    private int nrPts1; // nrPts1 x nrPts1

//    private static float s = 8f;
    private static int bg = 20; 	// background

//	private boolean randomizeLocations = false;

	private static float Deg2Rad = (float) (Math.PI/180);

	public Generator(int _nrPts1){nrPts1 = _nrPts1;}

	public ImagePlus runDisProportional(float snr, float diam1, float diam2, float diam3, float sc,
										File out_bif, File out_end, File out_non, File out_log, File out_img)
	{

        int nrPts2 = nrPts1*nrPts1;

        float diam_max = Math.max(diam1, Math.max(diam2, diam3));

        sc = (sc<=3)? 3 : sc ;

        float	branch_length = sc * diam_max;  // s*Dmax = branch_length
        float   branch_separation = branch_length;
        float	margin 	= branch_separation;    // pix
        float   separation_scale = sc/3;     // hardcoded
        int     minAngleDeg = (int) Math.round( 2f * Math.asin(1f/(2f*separation_scale)) * (180f/Math.PI) ); // Dmax branches are separable at separation_scale

        int W = Math.round(2 * margin + (2*nrPts1) * branch_separation + (nrPts1+1) * branch_separation);
        int H = W;

        // parameters used fo generate
        try {
            PrintWriter writer = new PrintWriter(out_log.getAbsolutePath());

            // first line of csv has header
            writer.print(String.format("SNR,\t p1,\t p2,\t p3,\t d1,\t d2,\t d3,\t s\n"));
            float pp1, pp2, pp3, pp;
            pp = diam1 + diam2 + diam3;
            pp1 = diam1 / pp;
            pp2 = diam2 / pp;
            pp3 = diam3 / pp;
            writer.print(String.format("%4.2f,\t %4.2f,\t %4.2f,\t %4.2f,\t %4.2f,\t %4.2f,\t %4.2f,\t %4.2f\n", snr, pp1, pp2, pp3, diam1, diam2, diam3, sc));
            writer.close();
        }
        catch (FileNotFoundException ex) {}

        PrintWriter writer_bif = null;
        PrintWriter writer_end = null;
        PrintWriter writer_non = null;

		// initialize file outputs
		try {
			writer_bif = new PrintWriter(new BufferedWriter(new FileWriter(out_bif, true)));
			writer_bif.println("#"); // dateFormat.format(date)
            writer_end = new PrintWriter(new BufferedWriter(new FileWriter(out_end, true)));
            writer_end.println("#"); // dateFormat.format(date)
			writer_non = new PrintWriter(new BufferedWriter(new FileWriter(out_non, true)));
			writer_non.println("#"); // dateFormat.format(date)
        } catch (IOException e) {}

		//GaussianBlur gauss 		= new GaussianBlur();
		PoissonGenerator poiss 	= new PoissonGenerator();  // noise

		float fg = foregroundLevel(bg, snr); // quadratic equation solving required snr

		Overlay ov = new Overlay(); // output image will have an overlay with gnd tth critpoints

		FloatProcessor outIp = new FloatProcessor(W, H);

		// initialize
		for (int j=0; j<W*H; j++) outIp.setf(j, 0);

        // form the skeleton of junction geometry
		float[][][] points = new float[nrPts2][4][2];
		int cnt = 0;
		for (int row = 0; row < nrPts1; row++) {
			for (int col = 0; col < nrPts1; col++) {
				float locx = margin + (col+1) * branch_separation + (2*col+1) * branch_length; //row*((W-2*margin)/(nrPts1-1));
				float locy = margin + (row+1) * branch_separation + (2*row+1) * branch_length; //col*((H-2*margin)/(nrPts1-1));
				points[cnt] = createBifPointsAt(locx, locy, branch_length, minAngleDeg); // 4(#points)x2(x,y)
				cnt++;
			}
		}

		for (int pIdx=0; pIdx<nrPts2; pIdx++) {  // loop all bifurcations

			// record ground truth in swc format
            int swc_idx = pIdx + 1;
            writer_bif.println(String.format("%d %1d %4.2f %4.2f %4.2f %4.2f %d", swc_idx, 5, points[pIdx][0][0], points[pIdx][0][1], 0f, diam_max, -1));

            writer_end.println(String.format("%d %1d %4.2f %4.2f %4.2f %4.2f %d", 3*swc_idx-2, 6, points[pIdx][1][0], points[pIdx][1][1], 0f, diam_max, -1));
            writer_end.println(String.format("%d %1d %4.2f %4.2f %4.2f %4.2f %d", 3*swc_idx-1, 6, points[pIdx][2][0], points[pIdx][2][1], 0f, diam_max, -1));
            writer_end.println(String.format("%d %1d %4.2f %4.2f %4.2f %4.2f %d", 3*swc_idx-0, 6, points[pIdx][3][0], points[pIdx][3][1], 0f, diam_max, -1));

            // add overlay
            OvalRoi pt = new OvalRoi(points[pIdx][0][0]+.5f-diam_max/2f, points[pIdx][0][1]+.5f-diam_max/2f, diam_max, diam_max);
            pt.setStrokeColor(Color.RED);
			ov.add(pt);

            for (int cnt_rest=1; cnt_rest<points[pIdx].length; cnt_rest++) {
                pt = new OvalRoi(points[pIdx][cnt_rest][0]+.5f-diam_max/2f, points[pIdx][cnt_rest][1]+.5f-diam_max/2f, diam_max, diam_max);
                pt.setStrokeColor(Color.YELLOW);
                pt.setFillColor(Color.YELLOW);
                ov.add(pt);
            }

			// set configuration
			float[] diameters = new float[3];
			diameters[0] = diam1;
			diameters[1] = diam2;
			diameters[2] = diam3;

			// draw junctions
			drawBif(points[pIdx], diameters, fg, outIp); // normalized [0, fg] gaussian profiles model intensity of the cross section

		}

		writer_bif.close();
		writer_end.close();
		IJ.log("exported: " + out_bif.getAbsolutePath());
		IJ.log("exported: " + out_end.getAbsolutePath());

		// PSF emulate (Gaussian blur on outSlice) - SKIPPED, OBJECTS ARE bigger than just spots
		//gauss.blurGaussian(outIp, 0.5, 0.5, 0.005);

		int count_non = 0;
		for (int j=0; j<W*H; j++) {

			// add background to be able to emulate poisson noise everywhere
			float currVal = outIp.getf(j);
			currVal += bg;
			// poisson
			outIp.setf(j, (float) poiss.next(currVal));
			//outIp.setf(j, currVal);

			// negatives are defined here, take the points that are far enough from critpoints and
			// take them with some probability

			if (Math.random()<=.05) {

				int loc_x = j - W * (j / W);
				int loc_y = j / W;

                if (loc_x>margin && loc_x<W-margin && loc_y>margin && loc_y<H-margin) {

				float min_d2 = Float.MAX_VALUE;

				for (int points_i=0; points_i<points.length; points_i++) {
					for (int points_k=0; points_k<points[points_i].length; points_k++) {

						float points_x = points[points_i][points_k][0];
						float points_y = points[points_i][points_k][1];

						float d2 = (float) (Math.pow(points_x-loc_x, 2)+Math.pow(points_y-loc_y, 2));

						if (d2 < min_d2) {
							min_d2 = d2;
						}

					}
				}

				if (min_d2>=Math.pow(2*diam_max, 2)){
					// far enough to add it as non point
					count_non++;
					writer_non.println(String.format("%d %1d %4.2f %4.2f %4.2f %4.2f %d", count_non, 7, (float)loc_x, (float)loc_y, 0f, 1f, -1));
					PointRoi pt = new PointRoi(loc_x+.5f, loc_y+.5f);
					pt.setFillColor(Color.BLUE);
					ov.add(pt);

				}

                }

			}

		}

		writer_non.close();
		IJ.log("exported: " + out_non.getAbsolutePath());

		ImagePlus outImp = new ImagePlus("SNR_"+IJ.d2s(snr,1)+",d1_"+IJ.d2s(diam1,1)+",d2_"+IJ.d2s(diam2,1)+",d3_"+IJ.d2s(diam3,1)+",s_"+IJ.d2s(sc,1), outIp);
		outImp.setOverlay(ov);

        // save it
        FileSaver fs = new FileSaver(outImp);
        fs.saveAsTiff(out_img.getAbsolutePath());
        IJ.log("exported: "+out_img.getAbsolutePath());
        IJ.log(nrPts2 + " bifs " + 4*nrPts2 + " ends " + count_non + " nons");

        return outImp;

	}

	private static float min(float[] in)
	{
		float out = in[0];
		for (int a=1; a<in.length; a++) {
			if (in[a]<out)
				out = in[a];
		}
		return out;
	}

	private static float max(float[] in)
	{
		float out = in[0];
		for (int a=1; a<in.length; a++) {
			if (in[a]>out)
				out = in[a];
		}
		return out;
	}

	private static float distToSegment(
												   int x, int y,   					// point considered
												   float[] pStart, float[] pEnd//,    // limits
												   //float diam                       //
	)
	{

		float d = 0;

		double[] p_b = new double[2];

		p_b[0] = x;
		p_b[1] = y;

		double[] n = new double[2];
		double nLen = Math.sqrt(Math.pow(pEnd[0]-pStart[0],2)+Math.pow(pEnd[1]-pStart[1],2));
		n[0] = (pEnd[0]-pStart[0])/nLen;
		n[1] = (pEnd[1]-pStart[1])/nLen;

		double proj = (p_b[0] - pEnd[0]) * (pStart[0]-pEnd[0]) + (p_b[1] - pEnd[1]) * (pStart[1]-pEnd[1]);
		if(proj<0){
			return (float) Math.sqrt(Math.pow(p_b[0]-pEnd[0],2)+Math.pow(p_b[1]-pEnd[1],2));
		}

		proj = (p_b[0] - pStart[0]) * n[0] + (p_b[1] - pStart[1]) * n[1];
		if(proj<0){
			return (float) Math.sqrt(Math.pow(p_b[0]-pStart[0],2)+Math.pow(p_b[1]-pStart[1],2));
		}

		p_b[0] = p_b[0] - pStart[0] - proj * n[0];
		p_b[1] = p_b[1] - pStart[1] - proj * n[1];

		return (float) Math.sqrt(p_b[0]*p_b[0]+p_b[1]*p_b[1]);


	}

	private static void drawBif(float[][] pts, float[] diams, float h, FloatProcessor image) // , float val
	{

		float extDiam = max(diams);

		int x1 = (int) Math.floor(
			Math.min(Math.min(pts[0][0], pts[1][0]), Math.min(pts[2][0], pts[3][0])) - extDiam
		);
		int x2 = (int) Math.ceil(
			Math.max(Math.max(pts[0][0], pts[1][0]), Math.max(pts[2][0], pts[3][0])) + extDiam
		);

		int y1 = (int) Math.floor(
			Math.min(Math.min(pts[0][1], pts[1][1]), Math.min(pts[2][1], pts[3][1])) - extDiam
		);
		int y2 = (int) Math.ceil(
			Math.max(Math.max(pts[0][1], pts[1][1]), Math.max(pts[2][1], pts[3][1])) + extDiam
		);

		x1 = (x1<0)? 0 : (x1>=image.getWidth())? image.getWidth()-1 : x1 ;
		x2 = (x2<0)? 0 : (x2>=image.getWidth())? image.getWidth()-1 : x2 ;

		y1 = (y1<0)? 0 : (y1>=image.getHeight())? image.getHeight()-1 : y1 ;
		y2 = (y2<0)? 0 : (y2>=image.getHeight())? image.getHeight()-1 : y2 ;

		for (int x=x1; x<=x2; x++) {
			for (int y=y1; y<=y2; y++) {

				for (int ptIdx = 1; ptIdx<=3; ptIdx++) {


					float dist = distToSegment(x, y, pts[0], pts[ptIdx]);

					// apply gaussian profile -> assume diam==2*std

//					float std = diams[ptIdx-1]/2; //*** HERE IT RELATES DIAMETER WITH STANDARD DEVIATION  DIAMETER ~ 2 * std
					float std = diams[ptIdx-1]/2.5f; //*** HERE IT RELATES DIAMETER WITH STANDARD DEVIATION  DIAMETER ~ 2.5 * std

					if (dist<=3*std) {
						float val = (float) (h * Math.exp(-(dist*dist)/(2*std*std))); // normalize height to h
						val = Math.max(val, image.getf(x, y));
						image.setf(x, y, val);
						//break;
					}
					//else {
						//float val = image.getf(x, y);
						//image.setf(x, y, val);
						//break;
					//}

				}

			}
		}


	}

	private static float foregroundLevel(float bg, float snr)
	{

//		return (float) (0.5*(2*bg+snr*snr+Math.sqrt(4*bg*snr+snr*snr)));
//		return (float) (0.5 * (snr*snr + snr * Math.sqrt(4*bg+1))); this was freaking wrong!!!


        // explanation b=background, a=added on top of b
        // a        signal
        // a+b      background
        // both poisson distribution (mean = lambda, variance = lambda, sigma = sqrt(lambda))
        // snr = mean(signal)/sigma(signal) = a / sqrt(a+b)
        // solve quadratic equation for a
        // quadratic equation is:
        // a^2  -snr^2  *a  -snr^2*b   = 0
        // a  = (1/2) (snr^2 + snr * sqrt(snr^2+4b))
        return (float) (0.5 * (snr * snr + snr * Math.sqrt(4 * bg + snr * snr)));

	}

	private static float[][] createRandomBifPoints(float minX, float maxX, float minY, float maxY, float branchLength)
	{
		float[][] p = new float[4][2];

		p[0][0] = random(minX, maxX); // col
		p[0][1] = random(minY, maxY); // row

		float refAng = random(0, 359);

		p[1][0] = (float) (p[0][0] +  branchLength * Math.cos(refAng * Deg2Rad)); // col
		p[1][1] = (float) (p[0][1] -  branchLength * Math.sin(refAng * Deg2Rad)); // row

		boolean confOK = false;
		int alfa1=0, alfa2=0, alfa3=0;

		while (!confOK) {

			alfa1 = Math.round(random(10, 180));//  ;
			alfa2 = Math.round(random(10, 180));//random(10, 180) * Deg2Rad ;
			alfa3 = 360-alfa1-alfa2;//random(10, 180) * Deg2Rad ;

			confOK = alfa3 >=10 && alfa3<=180;

		}

		p[2][0] = (float) (p[0][0] +  branchLength * Math.cos(Tools.wrap_360(refAng+alfa1) * Deg2Rad)); // col
		p[2][1] = (float) (p[0][1] -  branchLength * Math.sin(Tools.wrap_360(refAng+alfa1) * Deg2Rad)); // row

		p[3][0] = (float) (p[0][0] +  branchLength * Math.cos(Tools.wrap_360(refAng+alfa1+alfa2) * Deg2Rad));
		p[3][1] = (float) (p[0][1] -  branchLength * Math.sin(Tools.wrap_360(refAng+alfa1+alfa2) * Deg2Rad));

		return p;

	}

	private static float[][] createBifPointsAt(float atX, float atY, float branchLength, int minAngleDeg) // 4(#points)x2(x,y)
	{
		float[][] p = new float[4][2];

		p[0][0] = atX;//random(minX, maxX); // col
		p[0][1] = atY;//random(minY, maxY); // row

		float refAng = random(0, 359);

		p[1][0] = (float) (p[0][0] +  branchLength * Math.cos(refAng * Deg2Rad)); // col
		p[1][1] = (float) (p[0][1] -  branchLength * Math.sin(refAng * Deg2Rad)); // row

		boolean confOK = false;
		int alfa1=0, alfa2=0, alfa3=0;

		while (!confOK) {

			alfa1 = Math.round(random(minAngleDeg, 180));
			alfa2 = Math.round(random(minAngleDeg, 180));
			alfa3 = 360-alfa1-alfa2;//random(10, 180) * Deg2Rad ;

			confOK = alfa3 >=minAngleDeg && alfa3<=180;

		}

		p[2][0] = (float) (p[0][0] +  branchLength * Math.cos(Tools.wrap_360(refAng+alfa1) * Deg2Rad)); // col
		p[2][1] = (float) (p[0][1] -  branchLength * Math.sin(Tools.wrap_360(refAng+alfa1) * Deg2Rad)); // row

		p[3][0] = (float) (p[0][0] +  branchLength * Math.cos(Tools.wrap_360(refAng+alfa1+alfa2) * Deg2Rad));
		p[3][1] = (float) (p[0][1] -  branchLength * Math.sin(Tools.wrap_360(refAng+alfa1+alfa2) * Deg2Rad));

		return p;

	}

	private static float random(float min, float max)
	{
		Random rand = new Random();
		return min + rand.nextFloat()*(max-min);
	}

	private static int random(int min, int max)
	{
		Random rand = new Random();
		return min + rand.nextInt(max-min+1);
	}

	private float[] randomVec(float min, float max)
	{
		float[] vec = new float[3];
		vec[0] = random(min, max);
		vec[1] = random(min, max);
		vec[2] = random(min, max);
		return vec;
	}

	private float var3(float a, float b, float c)
	{
		float var = 0;
		float m = (a+b+c)/3;
		var = (1/3f) * ((a-m)*(a-m) + (b-m)*(b-m) + (c-m)*(c-m));
		return var;
	}

}