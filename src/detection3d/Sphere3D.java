package detection3d;

import detection.Interpolator;
import ij.ImageStack;
import ij.process.ByteProcessor;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 11/16/13
 * Time: 10:17 AM
 */
public class Sphere3D {

	// spherical coordinates
	// azimuth 		theta [0, 2PI)
	// polar angle  phi [PI/2(top), -PI/2(bottom)]

	private static float arcRes 	= 1.0f;
	private static float PI 		= (float) Math.PI;
	private static float TwoPI 		= (float) (2*Math.PI);
	private static float HalfPI 	= (float) (Math.PI/2);
	private static float samplingStep = 1.0f;
	public static float     R_FULL 	= 1.00f;
	public static float     T_HALF 	= 0.50f;
	public static int 		weightStdRatioToD = 4;  // could be a parameter

	private float 	radius;                         // sphere radius
	private static 	float 	arcNbhood = 1.5f;       // 2 * arcRes;
	private int 	W, H, N;                        // W,H are width and height of the 2d profile, N defines optimal sampling
	private int 	limR, limT;

	ArrayList<int[]>		vizXY = new ArrayList<int[]>();         // XY

	ArrayList<float[]> 		elems = new ArrayList<float[]>(); 	    // list of elements (phi, theta) covering the sphere

	ArrayList<int[]> 		masks = new ArrayList<int[]>(); 	    // list of list indexes of the neighbours for each list element

	ArrayList<float[][]> 	offstXYZ = new ArrayList<float[][]>(); 	// list of filter offsets for each direction

    float[] weights;

	public Sphere3D(float neuronDiam, float scale) {

		this.radius 	= scale * neuronDiam;
		this.N 	= (int) Math.ceil( ( (2 * Math.PI * radius) / arcRes - 1 ) / 4 );
		this.W  = 4*N + 1;
		this.H  = 2*N + 1;

		this.limR = (int) Math.ceil(R_FULL*neuronDiam/samplingStep);   // how many to take radially with given sampling step
		this.limT = (int) Math.ceil(T_HALF*neuronDiam/samplingStep);

		vizXY.clear();
		elems.clear();

		float stepPhi = PI / (H-1); // step when last included: []
		for (int phiIdx = 0; phiIdx<H; phiIdx++) {
			float currPhi = (PI/2) - phiIdx * stepPhi;

			int Wcurr = 2 * (int) (Math.round(Math.cos(currPhi) * 2 * N)) + 1;

			float stepTheta = TwoPI / Wcurr; // last not included: [)
			for (int thetaIdx = 0; thetaIdx<Wcurr; thetaIdx++) {

				float currTheta = thetaIdx * stepTheta;
				elems.add(new float[]{currPhi, currTheta});
				vizXY.add(new int[]{  (W-1)/2 - (Wcurr-1)/2 + thetaIdx  , phiIdx });

			}

		}

		System.out.println("profile size: "+elems.size()+" elements");

		// initialize mask list
		masks.clear();
		for (int ii=0; ii<elems.size(); ii++) {

			float phi = elems.get(ii)[0];
			float theta = elems.get(ii)[1];

			// loop the rest of the elements to see if they are in the angNeighbourhood range
			ArrayList<Integer> nbrs = new ArrayList<Integer>();
			for (int jj = 0; jj<elems.size(); jj++) {

				if (jj!=ii) {

					// check
					float phi_test = elems.get(jj)[0];
					float theta_test = elems.get(jj)[1];

					float ang_btw = arcBetweenDirections(phi, theta, phi_test, theta_test);
					if (ang_btw<=arcNbhood) {
						nbrs.add(jj);
					}

                }

			}

			// convert list to regular array and add
			int[] nbrsArray = new int[nbrs.size()];
			for (int iii=0; iii<nbrs.size(); iii++) {
				nbrsArray[iii] = nbrs.get(iii);
			}
			masks.add(nbrsArray);

		}

		// define coordinates for each direction
		offstXYZ.clear();
        // weights define only once
        weights = new float[(2*limT+1)*(2*limT+1)*(limR+1)];
        float sumWgt = 0;

		for (int ii = 0; ii<elems.size(); ii++) {

			float phi = elems.get(ii)[0];
			float theta	= elems.get(ii)[1];

			/*
				form sampling (offset) points
			 */

			float[][] offsetsPerDirection = new float[(2*limT+1)*(2*limT+1)*(limR+1)][3];

			int cnt = 0;

			for (int k=-limR; k<=0; k++) {

				for (int i=-limT; i<=limT; i++) {

					for (int j = -limT; j<=limT; j++) {

						float px = i * samplingStep;
						float py = j * samplingStep;
						float pz = k * samplingStep;

						offsetsPerDirection[cnt][0] = px;
						offsetsPerDirection[cnt][1] = py;
						offsetsPerDirection[cnt][2] = pz;

						//*** HERE IT DEFINES THE FILTER PROFILE WEIGHTS (only in first iteration, the rest are the same)
						if (ii==0) {
                            float dstAxis = point2line(0, 0, 0, 0, 0, -1, px, py, pz);
                            weights[cnt] = (float) Math.exp(-(dstAxis*dstAxis)/(2*(neuronDiam/weightStdRatioToD)*(neuronDiam/weightStdRatioToD)));
                            sumWgt += weights[cnt];
                        }

						cnt++;

					}

				}

			}

			/*
				transformation for offsets before adding
			 */
			transZ(radius, offsetsPerDirection);
			rotY(-phi+HalfPI, offsetsPerDirection);
			rotZ(theta, offsetsPerDirection);
			offstXYZ.add(offsetsPerDirection); //store

		}

        /*
				normalize weights
	    */
        for (int iii=0; iii<weights.length; iii++) {
            weights[iii] /= sumWgt;
        }
        // weights.add(weights);

	}

	private float angBetweenDirections(float phi1, float theta1, float phi2, float theta2) {

		// TODO: radius not necessary
		float x1 = getX(radius, phi1, theta1);
		float y1 = getY(radius, phi1, theta1);
		float z1 = getZ(radius, phi1);

		float x2 = getX(radius, phi2, theta2);
		float y2 = getY(radius, phi2, theta2);
		float z2 = getZ(radius, phi2);

		return (float) Math.acos( (x1*x2+y1*y2+z1*z2)/(radius * radius) );

	}

	private float arcBetweenDirections(float phi1, float theta1, float phi2, float theta2){

		float x1 = getX(radius, phi1, theta1);
		float y1 = getY(radius, phi1, theta1);
		float z1 = getZ(radius, phi1);

		float x2 = getX(radius, phi2, theta2);
		float y2 = getY(radius, phi2, theta2);
		float z2 = getZ(radius, phi2);

		return radius * (float) Math.acos( (x1*x2+y1*y2+z1*z2)/(radius * radius) );

	}

	private float point2line(float n1x, float n1y, float n1z,
							  float n2x, float n2y, float n2z,
							  float px, float py, float pz) {

		// line is defined with l and m

		float d = 0;

		double[] p_b = new double[3];

		// p - b
		//p_b[0] = px;// - b[0];
		//p_b[1] = py;// - b[1];
		//p_b[2] = pz;// - b[1];

		//double[] n21 = new double[3];
		float n21Len = (float) Math.sqrt(Math.pow(n2x-n1x,2)+Math.pow(n2y-n1y,2)+Math.pow(n2z-n1z,2));
		float n21x = (n2x-n1x)/n21Len;
		float n21y = (n2y-n1y)/n21Len;
		float n21z = (n2z-n1z)/n21Len;

		float proj = (px - n1x) * n21x + (py - n1y) * n21y + (pz - n1z) * n21z; // dot prod

//		if(Math.abs(proj)<Double.MIN_VALUE){
//			return Math.sqrt( Math.pow(p_b[0] - n2x, 2) + Math.pow(p_b[1] - n2y, 2) ); //Double.MAX_VALUE;
//		}
//
//		proj = (p_b[0] - n1x) * n[0] + (p_b[1] - n1y) * n[1];
//		if(Math.abs(proj)<Double.MIN_VALUE){
//			return Math.sqrt( Math.pow(p_b[0]-n1x, 2) + Math.pow(p_b[1] - n1y, 2)); //Double.MAX_VALUE;
//		}

		//IJ.log("nLen: "+nLen+" -> "+n[0]+","+n[1]+" proj: "+proj);

		p_b[0] = -(px - n1x) + proj * n21x;
		p_b[1] = -(py - n1y) + proj * n21y;
		p_b[2] = -(pz - n1z) + proj * n21z;

//		distance_from_line = vectorNorm(distance_2d);

		return (float) Math.sqrt(p_b[0]*p_b[0] + p_b[1]*p_b[1] + p_b[2]*p_b[2]);

	}

	private static float Rad2Deg(float angRad)
	{
		return (angRad/PI)*180f;
	}

	public ByteProcessor drawProfile(float[] profile) {

		// loops vizXY and corresponding profile values
		ByteProcessor outIp = new ByteProcessor(W, H);
		for (int ii=0; ii<elems.size(); ii++) {
			outIp.set(vizXY.get(ii)[0], vizXY.get(ii)[1], (byte) ((int)profile[ii] & 0xff));
		}

		return outIp;

	}

	public float[] extractProfile(float atX, float atY, float atZ, float[][][] img3d_zxy, float zDist) {
		// loops all elements and calculates filter outputs -> profile
		float[] out = new float[offstXYZ.size()];

		for (int profleIdx=0; profleIdx<offstXYZ.size(); profleIdx++) {

			float value = 0;

			for (int offsetIdx=0; offsetIdx<offstXYZ.get(profleIdx).length; offsetIdx++) {

				float x_offs_pix = offstXYZ.get(profleIdx)[offsetIdx][0];
				float y_offs_pix = offstXYZ.get(profleIdx)[offsetIdx][1];
				float z_offs_lay = offstXYZ.get(profleIdx)[offsetIdx][2] / zDist; // convert to layers

                float imgVal = Interpolator.interpolateAt(atX+x_offs_pix, atY+y_offs_pix, atZ+z_offs_lay, img3d_zxy);
				value += weights[offsetIdx] * imgVal;

			}

			out[profleIdx] = value;

		}

		return out;
	}

    private ArrayList<Integer> profilePeaks(float[] profile) {   // phi, theta

        // search for local maxima and corresponding (phi, theta)
        ArrayList<Integer> peaks = new ArrayList<Integer>();

        for (int direcIdx=0; direcIdx<profile.length; direcIdx++) {
            // calculate the index of local nbhood max for this profile sample

            boolean isLocMax = true;

            // check neighbours for this one
            for (int nbrLoop=0; nbrLoop<masks.get(direcIdx).length; nbrLoop++) {

                int nbrIdx = masks.get(direcIdx)[nbrLoop];
                if (profile[nbrIdx]>profile[direcIdx]) {
                    isLocMax = false;// check the profile value at the index
                }

            }

            if (isLocMax) {
                peaks.add(direcIdx); //new float[]{elems.get(direcIdx)[0], elems.get(direcIdx)[1]});
            }

        }

        return peaks;

    }

    public ArrayList<int[]> profilePeaksXY(float[] profile){ // profile peaks as planisphere image coordinates

        // indexes of peaks
        ArrayList<Integer> peakIdxs = profilePeaks(profile);
        // 2d coordinates of the peaks in planisphere representation
        ArrayList<int[]> out = new ArrayList<int[]>(peakIdxs.size());

        for (int t=0; t<peakIdxs.size(); t++) {
            out.add(vizXY.get(peakIdxs.get(t)));
        }

        return out;

    }

    public ArrayList<int[]> profilePeaksXYZ(float[] profile, float zDistImage) { // profile peaks as locations in 3d image (stacked set of 2d images)

        // indexes of peaks
        ArrayList<Integer> peakIdxs = profilePeaks(profile);

        // index can be automatically mapped to (x,y,z) coordinate

        // 3d coordinates (x[px],y[px],z[lay])
        ArrayList<int[]> out = new ArrayList<int[]>(peakIdxs.size());

        for (int t=0; t<peakIdxs.size(); t++) {

            // convert to xyz
            int indexOfPeakDirection = peakIdxs.get(t);

            int[] peaksXYZ = new int[3];

            float phi       = elems.get(indexOfPeakDirection)[0]; // position 0 iz marked as phi
            float theta     = elems.get(indexOfPeakDirection)[1]; // position 1 iz marked as phi

            // extract values in pixels
            peaksXYZ[0] = Math.round( getX(radius, phi, theta) );                 // x
            peaksXYZ[1] = Math.round( getY(radius, phi, theta) );                 // y
            peaksXYZ[2] = Math.round( getZ(radius, phi       ) / zDistImage );    // layer

            out.add(peaksXYZ);
        }

        return out;

    }

	private float getX(float r, float phi, float theta){return r * (float) Math.cos(phi) * (float) Math.cos(theta);}

	private float getY(float r, float phi, float theta){return r * (float) Math.cos(phi) * (float) Math.sin(theta);}

	private float getZ(float r, float phi             ){return r * (float) Math.sin(phi)                          ;}

	private void  transX(float tx, float[][] coords){
		for (int i=0; i<coords.length; i++) {
			coords[i][0] += tx;
		}
	}

	private void transY(float ty, float[][] coords){
		for (int i=0; i<coords.length; i++) {
			coords[i][1] += ty;
		}
	}

	private void transZ(float tz, float[][] coords){
		for (int i=0; i<coords.length; i++){
			coords[i][2] += tz;
		}
	}

	private void rotX(float ang, float[][] coords) {
		for (int i=0; i<coords.length; i++) {
			float y_temp = coords[i][1];
			float z_temp = coords[i][2];
			coords[i][1] = y_temp * (float) Math.cos(ang) - z_temp * (float) Math.sin(ang);
			coords[i][2] = y_temp * (float) Math.sin(ang) + z_temp * (float) Math.cos(ang);
		}

	}

	private void rotY(float ang, float[][] coords) {
		for (int i=0; i<coords.length; i++) {
			float x_temp = coords[i][0];
			float z_temp = coords[i][2];
			coords[i][0] = x_temp * (float) Math.cos(ang) + z_temp * (float) Math.sin(ang);
			coords[i][2] = -x_temp * (float) Math.sin(ang) + z_temp * (float) Math.cos(ang);
		}
	}

	private void rotZ(float ang, float[][] coords) {
		for (int i=0; i<coords.length; i++) {
			float x_temp = coords[i][0];
			float y_temp = coords[i][1];
			coords[i][0] = x_temp * (float) Math.cos(ang) - y_temp * (float) Math.sin(ang);
			coords[i][1] = x_temp * (float) Math.sin(ang) + y_temp * (float) Math.cos(ang);
		}
	}

	public void printMasks(String dirPath){

		String fileName = dirPath+File.separator+"S3D_nbhood_sampling_all.swc";

		PrintWriter logWriter = null;

		try {
			logWriter = new PrintWriter(fileName);
			logWriter.print("");
			logWriter.close();
		} catch (FileNotFoundException ex) {}

		try {
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(fileName, true)));
			logWriter.println("# all nbhoods");
		} catch (IOException e) {}

		int cnt = 1;
		// take every third (just to check & visualize)
		for (int elemsIdx=0; elemsIdx<elems.size(); elemsIdx+=20) {

			float phi 	= elems.get(elemsIdx)[0];
			float theta	= elems.get(elemsIdx)[1];

			// export swc file with sampling for this direction only
			String fileName1 = dirPath+File.separator+"S3D_nbhood_sampling_"+elemsIdx+".swc";

			PrintWriter logWriter1 = null;

			try {
				logWriter1 = new PrintWriter(fileName1);
				logWriter1.print("");
				logWriter1.close();
			} catch (FileNotFoundException ex) {}

			try {
				logWriter1 = new PrintWriter(new BufferedWriter(new FileWriter(fileName1, true)));
				logWriter1.println("# filter "+elemsIdx);
			} catch (IOException e) {}

			int cnt1 = 1;

			for (int ii=0; ii<masks.get(elemsIdx).length; ii++){

				int maskElementIdx = masks.get(elemsIdx)[ii];
				float maskElementPhi = elems.get(maskElementIdx)[0];
				float maskElementTheta = elems.get(maskElementIdx)[1];

				// all filters log
				logWriter.println(String.format("%d 2 %4.2f %4.2f %4.2f 0.1 -1",
													   cnt++,
													   getX(radius, maskElementPhi, maskElementTheta),
													   getY(radius, maskElementPhi, maskElementTheta),
													   getZ(radius, maskElementPhi)));
				// individual filter log
				logWriter1.println(String.format("%d 2 %4.2f %4.2f %4.2f 0.1 -1",
														cnt1++,
														getX(radius, maskElementPhi, maskElementTheta),
														getY(radius, maskElementPhi, maskElementTheta),
														getZ(radius, maskElementPhi)));
			}

			// all nbhoods center
			logWriter.println(String.format("%d 3 %4.2f %4.2f %4.2f 0.50 -1",
												   cnt++,
												   0f,
												   0f,
												   0f));

			// individual nbhood center
			logWriter1.println(String.format("%d 3 %4.2f %4.2f %4.2f 0.50 -1",
													cnt1++,
													0f,
													0f,
													0f));

			// all nbhoods endpoint
			logWriter.println(String.format("%d 4 %4.2f %4.2f %4.2f 0.25 -1",
												   cnt++,
												   getX(radius, phi, theta),
												   getY(radius, phi, theta),
												   getZ(radius, phi)));

			// individual nbhood endpoint
			logWriter1.println(String.format("%d 4 %4.2f %4.2f %4.2f 0.25 -1",
													cnt1++,
													getX(radius, phi, theta),
													getY(radius, phi, theta),
													getZ(radius, phi)));

			logWriter1.close();
			System.out.println(fileName1+ "\t saved.");



		}

		logWriter.close();
		System.out.println(fileName+ "\t saved.");

	}

    public void printFilterDirections(String fileName){

        PrintWriter logWriter = null;

        try {
            logWriter = new PrintWriter(fileName);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}

        try {
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(fileName, true)));
            logWriter.println("# filter directions");
        } catch (IOException e) {}

        for (int ii=0; ii<elems.size(); ii++) {

            float phi 	= elems.get(ii)[0];
            float theta = elems.get(ii)[1];

            logWriter.println(String.format("%d 2 %4.2f %4.2f %4.2f 0.1 -1",
                    (ii + 1),
                    getX(radius, phi, theta),
                    getY(radius, phi, theta),
                    getZ(radius, phi)));

        }

        logWriter.println(String.format("%d 1 %4.2f %4.2f %4.2f 1.0 -1", (elems.size() + 1), 0f, 0f, 0f));
        logWriter.close();

    }


    public void printFilterOffsets(String dirPath){

		String fileName = dirPath+File.separator+"S3D_filter_sampling_all.swc";

		PrintWriter logWriter = null;

		try {
			logWriter = new PrintWriter(fileName);
			logWriter.print("");
			logWriter.close();
		} catch (FileNotFoundException ex) {}

		try {
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(fileName, true)));
			logWriter.println("# all filters");
		} catch (IOException e) {}

		int cnt = 1;
		// take every third (just to check & visualize)
		for (int elemsIdx=0; elemsIdx<elems.size(); elemsIdx+=20) {

			float phi 	= elems.get(elemsIdx)[0];
			float theta	= elems.get(elemsIdx)[1];

			// export swc file with sampling for this direction only
			String fileName1 = dirPath+File.separator+"S3D_filter_sampling_"+elemsIdx+".swc";

			PrintWriter logWriter1 = null;

			try {
				logWriter1 = new PrintWriter(fileName1);
				logWriter1.print("");
				logWriter1.close();
			} catch (FileNotFoundException ex) {}

			try {
				logWriter1 = new PrintWriter(new BufferedWriter(new FileWriter(fileName1, true)));
				logWriter1.println("# filter "+elemsIdx);
			} catch (IOException e) {}

			int cnt1 = 1;


			for (int ii=0; ii<offstXYZ.get(elemsIdx).length; ii++){

				// all filters log
				logWriter.println(String.format("%d 2 %4.2f %4.2f %4.2f 0.1 -1",
													   cnt++,
													   offstXYZ.get(elemsIdx)[ii][0],
													   offstXYZ.get(elemsIdx)[ii][1],
													   offstXYZ.get(elemsIdx)[ii][2]));
				// individual filter log
				logWriter1.println(String.format("%d 2 %4.2f %4.2f %4.2f 0.1 -1",
													   cnt1++,
													   offstXYZ.get(elemsIdx)[ii][0],
													   offstXYZ.get(elemsIdx)[ii][1],
													   offstXYZ.get(elemsIdx)[ii][2]));
			}

			// all filters center
			logWriter.println(String.format("%d 3 %4.2f %4.2f %4.2f 0.50 -1",
												   cnt++,
												   0f,
												   0f,
												   0f));

			// individual filter center
			logWriter1.println(String.format("%d 3 %4.2f %4.2f %4.2f 0.50 -1",
												   cnt1++,
												   0f,
												   0f,
												   0f));

			// all filters endpoint
			logWriter.println(String.format("%d 4 %4.2f %4.2f %4.2f 0.25 -1",
													cnt++,
													getX(radius, phi, theta),
													getY(radius, phi, theta),
													getZ(radius, phi)));

			// individual filter endpoint
			logWriter1.println(String.format("%d 4 %4.2f %4.2f %4.2f 0.25 -1",
												   cnt1++,
												   getX(radius, phi, theta),
												   getY(radius, phi, theta),
												   getZ(radius, phi)));

			logWriter1.close();
			System.out.println(fileName1+ "\t saved.");


		}

		logWriter.close();
		System.out.println(fileName+ "\t saved.");

	}

    public ImageStack visualizeMasks() {

        // stacked profiles with neighbourhood masks
        ImageStack isMasks = new ImageStack(W, H);

        // dummy profile with all values set to 255
        float[] profile = new float[W*H];
        Arrays.fill(profile, 100f);

        // every 3d direction will have its layer
        for (int ii=0; ii<masks.size(); ii++) {

            ByteProcessor layerMask = drawProfile(profile); // flat profile

            // central location
            int centerXplot = vizXY.get(ii)[0];
            int centerYplot = vizXY.get(ii)[1];
            layerMask.setf(centerXplot, centerYplot, 255f);

            // fill the location and its neighbours with 255
            for (int jj = 0; jj<masks.get(ii).length; jj++) {
                int neighbourIdx = masks.get(ii)[jj];
                int neighbourXplot = vizXY.get(neighbourIdx)[0];
                int neighbourYplot = vizXY.get(neighbourIdx)[1];
                layerMask.setf(neighbourXplot, neighbourYplot, 255f);
            }

            isMasks.addSlice(layerMask);

        }

        return isMasks;

    }

}
