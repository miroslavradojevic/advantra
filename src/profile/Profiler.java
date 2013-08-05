package profile;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/7/13
 * Time: 5:29 PM
 */
public class Profiler extends Thread {

	private int begN, endN;

	public static double 	neuronDiam;
	public static double 	scale;
	public static int 		outerRange;

	public static int 		resolDeg;

	public static ArrayList<ArrayList<double[]>> offsets; 	// each offset location score calculation will be parallelled
	public static ArrayList<ArrayList<Double>> weights; 	// weight for each offset - to emphasize body direction

	public static FloatProcessor    inip;
	public static ByteProcessor     maskip;

	public static float[][]     profiles;     // array of profiles
    public static int[][]       locations;    // [x, y]

	public static double    samplingStep = 0.5;
	public static float     TwoPI = (float) (Math.PI*2);

	public Profiler (int n0, int n1) {
		this.begN = n0;
		this.endN = n1;
	}

    public static void loadTemplate(ImageProcessor inip1, ByteProcessor maskip1) {

        /*
        set inip
         */
        inip = new FloatProcessor(inip1.getWidth(), inip1.getHeight());
        for (int i=0; i<inip1.getWidth()*inip1.getHeight(); i++) {
            inip.setf(i, inip1.getf(i));
        }

        /*
        set maskip
         */
        int cnt = 0;
        maskip = new ByteProcessor(maskip1.getWidth(), maskip1.getHeight());
        byte[] arr = (byte[]) maskip1.getPixels();
        for (int i=0; i<arr.length; i++) {
            if (arr[i]==(byte)255) {
                maskip.set(i, 255);
                cnt++;
            }
            else {
                maskip.set(i, 0);
            }
        }

        /*
        set locations
         */
        locations = new int[cnt][2];
        cnt = 0;
        for (int i=0; i<arr.length; i++) {
            if (arr[i]==(byte)255) {
                locations[cnt][0] = i%maskip.getWidth();
                locations[cnt][1] = i/maskip.getWidth();
                cnt++;
            }
        }

    }

	public static void loadParams(double neuronDiam1, double scale1, boolean showSampling) {

		neuronDiam 	= neuronDiam1;
		scale 		= scale1;
		outerRange 	= (int)Math.round(Math.sqrt(Math.pow(neuronDiam*scale+neuronDiam, 2) + Math.pow(neuronDiam, 2)));

        /*
        set offsets
         */
        double[] 	n1, n2, n0;
		n1 = new double[2];
		n2 = new double[2];
		n0 = new double[2];

		resolDeg = (int) ( Math.round( ( 2*Math.atan(1f / (2 * scale)) * (1f/8) / TwoPI) * 360 ) );
		resolDeg = (resolDeg>=1)? resolDeg : 1;

        // define width & length of the rectangular profile, width = 2*neuronDiam, height = 2*neuronDiam
		double r1 = neuronDiam * scale-neuronDiam;
		double r2 = neuronDiam * scale+neuronDiam;
		double r0 = neuronDiam * scale;

		/*
		create offsets, weights
		 */
		offsets = new ArrayList<ArrayList<double[]>>();
		weights = new ArrayList<ArrayList<Double>>();

		// +/-limR, +/-limT (used to limit index)
		int limR = (int) Math.ceil(neuronDiam/samplingStep);
		int limT = (int) Math.ceil(0.5*neuronDiam/samplingStep);

		/* 	these will be used to visualize the sampling and the shape of the profiles
			if showSampling was set to true
		*/
		Overlay ov  = new Overlay();
		ImageStack stackSampling = new ImageStack();
		ImageStack stackProfile  = new ImageStack();
		FloatProcessor stackProfileSlice = new FloatProcessor(2*limT+1, 2*limR+1);
		PointRoi pt;

		if (showSampling) {
			stackSampling 	= new ImageStack(2*outerRange+1, 2*outerRange+1);
			stackProfile 	= new ImageStack(2*limT+1, 2*limR+1);
			pt = new PointRoi(outerRange+0.5, outerRange+0.5);
			ov.add(pt);
		}

		for (int aDeg = 0; aDeg<360; aDeg+=resolDeg) {

			float aRad = ((float)aDeg/360)*TwoPI;

			double sumWgt = 0;
			ArrayList<double[]> offsetsAngleLoc = new ArrayList<double[]>();
			ArrayList<Double> 	offsetsAngleWgt = new ArrayList<Double>();

			n1[0] = r1 * (float) Math.cos(aRad);
			n1[1] = r1 * (-1) * (float) Math.sin(aRad);

			n2[0] = r2 * (float) Math.cos(aRad);
			n2[1] = r2 * (-1) * (float) Math.sin(aRad);

			n0[0] = r0 * (float) Math.cos(aRad);
			n0[1] = r0 * (-1) * (float) Math.sin(aRad);

			double dx = samplingStep*Math.sin(aRad);
			double dy = samplingStep*Math.cos(aRad);

			if (showSampling) {
				stackProfileSlice = new FloatProcessor(2*limT+1, 2*limR+1);
			}

			for (int i=-limR; i<=limR; i++) {

				for (int j=-limT; j<=limT; j++) {

					double px = n0[0] + j * dx + i * (-dy);
					double py = n0[1] + j * dy + i * dx;

					double dst = point2line(n1[0], n1[1], n2[0], n2[1], px, py);
					offsetsAngleLoc.add(new double[]{px, py});
					double weight = Math.exp(-(dst*dst)/(2*(neuronDiam/2)*(neuronDiam/2)));
					offsetsAngleWgt.add(weight);
					sumWgt += weight;

					if (showSampling) {

						pt = new PointRoi(outerRange+px+0.5, outerRange+py+0.5);
						pt.setPosition(stackSampling.getSize()+1);
						pt.setStrokeColor(Color.ORANGE);
						ov.add(pt);

						stackProfileSlice.setf(j+limT, i+limR, (float) weight);

					}

				}
			}

			if (showSampling) {
				stackSampling.addSlice(new ByteProcessor(2 * outerRange + 1, 2 * outerRange + 1));
				stackProfile.addSlice(stackProfileSlice);
			}

			// normalize
			for (int k=0; k<offsetsAngleWgt.size(); k++) {
				double newVal = offsetsAngleWgt.get(k) / sumWgt;
				offsetsAngleWgt.set(k, newVal);
			}

			offsets.add(offsetsAngleLoc);
			weights.add(offsetsAngleWgt);
		}

		if (showSampling) {

			ImagePlus sampling = new ImagePlus("sampling_scheme", stackSampling);
			sampling.setOverlay(ov);
			sampling.show();

			ImagePlus prof = new ImagePlus("filter_profiles", stackProfile);
			prof.show();

		}

		/*
		 allocate profiles container
		*/
		profiles = new float [locations.length][offsets.size()];

	}

	public static float[] extractProfile(int atX, int atY) { // profile at one location

		float[] profileOut = new float[offsets.size()];

		// calculate profile
		for (int offsetIdx = 0; offsetIdx < offsets.size(); offsetIdx++) {

			profileOut[offsetIdx] = 0;

			// calculate weighted response
			for (int k=0; k<offsets.get(offsetIdx).size(); k++) {
				profileOut[offsetIdx] += Interpolator.interpolateAt(atX+offsets.get(offsetIdx).get(k)[0],
																	atY+offsets.get(offsetIdx).get(k)[1],
																	inip
				) * weights.get(offsetIdx).get(k);
			}

			//profileOut[offsetIdx] /= offsets.get(offsetIdx).size();

		}

		return profileOut;
	}

    public void run(){ // considers begN and endN

        byte[] arr = (byte[]) maskip.getPixels();
        int locIdxProfile = 0;

        for (int locIdx=0; locIdx<arr.length; locIdx++) {
            if (arr[locIdx]==(byte)255) {

                int atX = locIdx%maskip.getWidth();
                int atY = locIdx/maskip.getWidth();

                for (int offsetIdx = begN; offsetIdx < endN; offsetIdx++) {

                    profiles[locIdxProfile][offsetIdx] = 0;
                    // calculate average for locIdx taking values from offsets(offsetIdx)
                    for (int k=0; k<offsets.get(offsetIdx).size(); k++) {
                        profiles[locIdxProfile][offsetIdx] += Interpolator.interpolateAt(
                                					atX+offsets.get(offsetIdx).get(k)[0],
                                					atY+offsets.get(offsetIdx).get(k)[1],
                                					inip
                        ) * weights.get(offsetIdx).get(k);
                    }

                    //profiles[locIdxProfile][offsetIdx] /= offsets.get(offsetIdx).size();

                }

                locIdxProfile++;   // increment if it was in the mask

            }
        }

    }

	private static double 	point2line(
										double n1x, 		// limit
										double n1y,

										double n2x, 		// limit
										double n2y,

										double px,        // point considered
										double py
	)
	{
		// line is defined with n1 and n2

		float d = 0;

		double[] p_b = new double[2];

		// p - b
		p_b[0] = px;// - b[0];
		p_b[1] = py;// - b[1];

		double[] n = new double[2];
		double nLen = Math.sqrt(Math.pow(n2x-n1x,2)+Math.pow(n2y-n1y,2));
		n[0] = (n2x-n1x)/nLen;
		n[1] = (n2y-n1y)/nLen;

		double proj = (p_b[0] - n2x) * (n1x-n2x) + (p_b[1] - n2y) * (n1y-n2y); // dot prod

		if(Math.abs(proj)<Double.MIN_VALUE){
			return Math.sqrt( Math.pow(p_b[0] - n2x, 2) + Math.pow(p_b[1] - n2y, 2) ); //Double.MAX_VALUE;
		}

		proj = (p_b[0] - n1x) * n[0] + (p_b[1] - n1y) * n[1];
		if(Math.abs(proj)<Double.MIN_VALUE){
			return Math.sqrt( Math.pow(p_b[0]-n1x, 2) + Math.pow(p_b[1] - n1y, 2)); //Double.MAX_VALUE;
		}

		//IJ.log("nLen: "+nLen+" -> "+n[0]+","+n[1]+" proj: "+proj);

		p_b[0] = p_b[0] - n1x - proj * n[0];
		p_b[1] = p_b[1] - n1y - proj * n[1];

//		distance_from_line = vectorNorm(distance_2d);

		return Math.sqrt(p_b[0]*p_b[0]+p_b[1]*p_b[1]);
	}

}