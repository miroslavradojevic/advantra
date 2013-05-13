package advantra.feature;

import ij.IJ;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class CircularConfiguration {

    public int			nrPeaks;

	public int[] 		angResDeg;
	public double[]		angResRad;

	public int[] 		angBtwPeakDeg;
	public double[] 	angBtwPeakRad;
	
	//public int 			nrRot;
	//public double 		rotStepRad;
	
//	public double 		minAngResRad; // only used to define how many rotations there are
//    public int          minAngResDeg; // only used to define how many rotations there are

	public double[] 	peaksRad;
	public double		r;

    public static float TwoPi = (float) (2*Math.PI);

	double[] 			filt;
	
	public CircularConfiguration(
            int[] angular_resolution_deg,
            int[] angles_between_peaks_deg,
			double radius
    )
	{

		r = radius;

		// filter used when calculating score(), only to initialize
		filt = new double[1];    // length will depend on th profile

		angResDeg = new int[angular_resolution_deg.length];
        for (int i = 0; i < angular_resolution_deg.length; i++){
            angResDeg[i] = angular_resolution_deg[i];
        }

		angResRad = new double[angular_resolution_deg.length];
        for (int i = 0; i < angular_resolution_deg.length; i++){
            angResRad[i] = (angResDeg[i]/360f)*TwoPi;
        }

		angBtwPeakDeg = new int[angles_between_peaks_deg.length];
		for (int i = 0; i < angles_between_peaks_deg.length; i++) {
			angBtwPeakDeg[i] = angles_between_peaks_deg[i];
		}
		
		angBtwPeakRad = new double[angles_between_peaks_deg.length];
		for (int i = 0; i < angles_between_peaks_deg.length; i++) {
			angBtwPeakRad[i] = (angles_between_peaks_deg[i]/360.0)*TwoPi;
		}

		nrPeaks = angles_between_peaks_deg.length;

        // find the smallest angular resolution to guide rotations
//        minAngResDeg = 360;
//        for (int i = 0; i < angResDeg.length; i++){
//            if (angResDeg[i]<minAngResDeg){
//                minAngResDeg = angResDeg[i];
//            }
//        }
//        minAngResRad = (minAngResDeg/360f)*TwoPi;
//		nrRot = (int)Math.round( TwoPi / (minAngResRad/2) );
//		rotStepRad = minAngResRad/2;

		peaksRad = new double[nrPeaks];
		// each configuration will be determined with
        // - rotations,
        // - peaks (and their angular width),
        // - angles between peaks
        // they will influence filter design
        // for this configuration
        // this is the core information of the config
        // where the peaks are
		
		//for (int cnt_rots = 0; cnt_rots < nrRot; cnt_rots++) {
			
//			double start_pos = 0;//cnt_rots*rotStepRad;

		peaksRad[0] 	= 0;

			for (int cnt_pks = 1; cnt_pks < nrPeaks; cnt_pks++) {
				
//				if(cnt_pks==0)
//					peaksRad[cnt_pks] 	= 0;//start_pos;
//				else
					peaksRad[cnt_pks] 	=
                            peaksRad[cnt_pks-1] +
                            angBtwPeakRad[cnt_pks-1];
				
			}
			
		//}

	}

	public void initFilter(int length)
	{
		filt = new double[length];
	}

    public void print()
    {

        System.out.print("circ. conf.  ");
        for (int j = 0; j < nrPeaks; j++){
            System.out.print(
                            IJ.d2s((peaksRad[j] / TwoPi) * 360, 2)+
                            " ( w "+
                            IJ.d2s((angResRad[j]/TwoPi)*360 ,2)+
                            " ) ; r "+
							IJ.d2s(r,2)
            );
        }

    }

    public float score(
            float[] val,
            float[] ang,
			float[] rad
    )
	{

        /*
        calculate score on configuration as
        highest score on all rotated versions of
        that configuration
         */

//        float 	score = Float.NEGATIVE_INFINITY;
		int sumPos, sumNeg;

        //for (int r = 0; r < nrRot; r++) {

            // create filt[]
            sumPos = sumNeg = 0;

            for (int i = 0; i < val.length; i++){

				if (rad[i]<=r){
					if(isOn(ang[i])){
						filt[i] = +1;
						sumPos++;
					}
					else{
						filt[i] = -1;
						sumNeg++;
					}
				}
				else{
					filt[i] = 0;
				}

            }

            // normalize filt[]
            for (int i1 = 0; i1 < filt.length; i1++) {
                if (filt[i1]>0) {
                    filt[i1] = filt[i1] / sumPos;
                }
				else if (filt[i1]<0){
					filt[i1] = filt[i1] / sumNeg;
				}
            }

            // score for this filt[]
            float sc = 0;
            for (int i = 0; i < filt.length; i++) {
                sc += filt[i]*val[i];
            }

//            if(sc>score || r==0){
//                score = sc;
//            }

        //}

        return sc;

    }
	
	private boolean isOn(float angle)
	{

			for (int p = 0; p < nrPeaks; p++) {
				if(Math.abs(wrap_PI(angle-peaksRad[p])) <= angResRad[p]/2){
					return true;
				}
			}

		return false;

	}
	
	public ImageProcessor plot()
	{

        int N = 101;
//		ImageStack viz_is = new ImageStack(N, N);

        int centerX = (N-1)/2;
        int centerY = (N-1)/2;
		int R2		= (N-1)/2;
		int R1      = (int) (r*R2);
//		for (int rotIdx = 0; rotIdx < nrRot; rotIdx++) {

	    ImageProcessor fp = new FloatProcessor(N, N);

            // fill the values
            for (int c = 0; c < fp.getWidth(); c++){
                for (int r = 0; r < fp.getHeight(); r++){
                    float ro, theta;
                    ro = (int)Math.sqrt(Math.pow((c-centerX), 2)+Math.pow((r-centerY), 2));

                    if (ro<=R1){

                        theta = (float) (Math.atan2(r-centerY, c-centerX) + Math.PI);

						if(isOn(theta)){

							fp.setf(c, r, +1);

						}
						else {

							fp.setf(c, r, -1);

						}

                    }

                }
            }

//            viz_is.addSlice("", fp);
//		}
		
		return fp;

	}

	public int[] getConfiguration()
	{
		int[] out_conf = new int[2*nrPeaks];
		for (int i = 0; i < nrPeaks; i++){
			out_conf[2*i] = angResDeg[i];
			out_conf[2*i+1] = angBtwPeakDeg[i];
		}
		return out_conf;
	}

	private static double wrap_PI(double in)
	{
		
		double out = in;
		
		while(out<=-Math.PI){
			out += 2*Math.PI;
		}
		while(out>Math.PI){
			out -= 2*Math.PI;
		}
		
		return out;
	}
	
}



//	public ImageStack plotFilter()
//	{
//
//		ImageStack viz_is = new ImageStack(400, 200);
//
//		int angular_res = 128;
//		int radius_res 	= 128;
//
//		float[] angles = new float[angular_res];
//		float[] radiuses = new float[radius_res];
//
//		for (int i = 0; i < angular_res; i++) {
//			angles[i] = i * (TwoPi / angular_res);
//		}
//
//		for (int i = 0; i < radius_res; i++) {
//			radiuses[i] = i * (1f / radius_res);
//		}
//
//
//////		for (int rotIdx = 0; rotIdx < nrRot; rotIdx++) {
////
////			float[] filt = new float[res];
////			float[] angles = new float[res];
////
////			int sumPos = 0;
////			int sumNeg = 0;
////
////			for (int i = 0; i < res; i++) {
////
////				angles[i] = i * (TwoPi / res);
////
////				boolean isON = isOn(angles[i]);
//////				// check if it belongs to ON
//////				for (int p = 0; p < nrPeaks; p++) {
//////					if(Math.abs(wrap_PI(angle-peaksRad[rotIdx][p])) <= angResRad[p]/2){
//////						isON = true;
//////						break;
//////					}
//////				}
////
////				if(isON){
////					filt[i] = +1;
////					sumPos++;
////				}
////				else{
////					filt[i] = -1;
////					sumNeg++;
////				}
////
////			}
////
////			// normalize filt[]
////			for (int i = 0; i < filt.length; i++) {
////				if(filt[i]>0){
////					filt[i] = filt[i] / (float)sumPos;// * ((float)sumNeg/(float)sumPos);
////				}
////                else{
////                    filt[i] = filt[i] / (float)sumNeg;
////                }
////			}
//
//			Plot p = new Plot("filter","angle","filter_value");
//		    p.addPoints(, filt, Plot.BOX);
//
//			p.setSize(400, 200);
//
//			viz_is.addSlice("rot."+0+"/"+nrRot, p.getProcessor());
//
////		}
//
//		return viz_is;
//
//	}


//	public double calculateScore(double[] profile_2PI){
//
//		/* calculate score on the configuration as
//		 * highest score on all rotated versions of
//		 * that configuration to be rotationally invariant
//		 */
//
//		double 	score = Double.NEGATIVE_INFINITY;
//		double[] filt = new double[profile_2PI.length];
//
//		int sumPos, sumNeg;
//
//		for (int r = 0; r < nrRot; r++) {
//
//			// create filt[]
//			sumPos = sumNeg = 0;
//
//			for (int i = 0; i < profile_2PI.length; i++) {
//
//				double angle = i*((2*Math.PI)/profile_2PI.length);
//
//                // check if it belongs to ON
//				boolean isON = false;
//				for (int p = 0; p < nrPeaks; p++) {
//					if(Math.abs(wrap_PI(angle-peaksRad[r][p])) <= wdthRad/2){
//						isON = true;
//						break;
//					}
//				}
//
//				if(isON){
//					filt[i] = +1;
//					sumPos++;
//				}
//				else{
//					filt[i] = -1;
//					sumNeg++;
//				}
//
//			}
//
//			// normalize filt[]
//			for (int i = 0; i < filt.length; i++) {
//				if(filt[i]>0){
//					filt[i] = filt[i] * ((float)sumNeg/(float)sumPos);
//				}
//			}
//
//			// score for this filt[]
//			double sc = 0;
//			for (int i = 0; i < filt.length; i++) {
//				sc += filt[i]*profile_2PI[i];
//			}
//
//			if(sc>score || r==0){
//				score = sc;
//			}
//
//		}
//
////		double[] xValues = new double[filt.length];
////		for (int i = 0; i < filt.length; i++) {
////			xValues[i] = i*(360f/filt.length);
////		}
//
////		Plot p = new Plot("", "ang [deg]", "");
////		p.setLimits(0, 360, -1.1, max_filt);
////		p.addPoints(xValues, profile_2PI, Plot.X);
////		p.addPoints(xValues, best_filt, Plot.LINE);
////		p.setSize(400, 200);
////		p.show();
//
//		return score;
//	}

//	public double calculateScore(float[] profile_2PI){
//
//		double 	score 	= Double.NEGATIVE_INFINITY;
//		float[] filt 	= new float[profile_2PI.length];
//
//		int sumPos, sumNeg;
//
//		for (int r = 0; r < nrRot; r++) {
//
//			sumPos = sumNeg = 0;
//
//			for (int i = 0; i < profile_2PI.length; i++) {
//
//				double angle = i*((2*Math.PI)/profile_2PI.length);
//
//				boolean isON = false;
//				for (int p = 0; p < nrPeaks; p++) {
//					if(Math.abs(wrap_PI(angle-peaksRad[r][p])) <= wdthRad/2){
//						isON = true;
//						break;
//					}
//				}
//
//				if(isON){
//					filt[i] = +1;
//					sumPos++;
//				}
//				else{
//					filt[i] = -1;
//					sumNeg++;
//				}
//
//			}
//
//			// normalize filt[]
//			for (int i = 0; i < filt.length; i++) {
//				if(filt[i]>0){
//					filt[i] = filt[i] * ((float)sumNeg/(float)sumPos);
//				}
//			}
//
//			// score for this filt[]
//			double sc = 0;
//			for (int i = 0; i < filt.length; i++) {
//				sc += filt[i]*profile_2PI[i];
//			}
//
//			if(sc>score || r==0){
//				score = sc;
//			}
//
//		}
//
//		return score;
//	}
//			float[] x = new float[plotResolution];
//			float[] y = new float[plotResolution];
//
//			for (int k = 0; k < plotResolution; k++) {
//
//				double thetaRad 	= k*((2*Math.PI)/plotResolution);
//				double ro = 1.5; // OFF
//
//				for (int l = 0; l < nrPeaks; l++) {
//					if(Math.abs(wrap_PI(thetaRad-peaksRad[rotIdx][l])) <= angResRad[l]/2){
//						ro=2.5;
//						break;
//					}
//				}
//
//				x[k] 	= (float) (ro * Math.sin(thetaRad));
//				y[k]	= (float) (ro * Math.cos(thetaRad));
//
//			}
//			int A = (int)(360.0/(float)minAngResDeg);
//			double[] x1 = new double[A];
//			double[] y1 = new double[A];
//			for (int k = 0; k < A; k++) {
//				x1[k] = 2 * Math.sin(k*minAngResDeg);
//				y1[k] = 2 * Math.cos(k*minAngResDeg);
//			}

//			Plot p = new Plot("profile","x","y", x, y);
//			p.setSize(400, 400);
//			p.setLimits(-2.5, 2.5, -2.5, 2.5);
//			p.setJustification(Plot.CENTER);
////			p.addPoints(x1, y1, Plot.BOX);