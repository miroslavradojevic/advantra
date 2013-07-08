package profile;

import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

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
	public static double 	outerRange;

	public static int 		resolDeg;

	public static ArrayList<ArrayList<double[]>> offsets;

	public static FloatProcessor    inip;
	public static ByteProcessor     maskip;

	public static float[][]     profiles;
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

	public static void loadParams(double neuronDiam1, double scale1) {

		neuronDiam = neuronDiam1;
		scale = scale1;
		outerRange = Math.sqrt(Math.pow(neuronDiam*scale+neuronDiam/2, 2) + Math.pow(neuronDiam/2, 2));

        /*
        set offsets
         */
        double[] 	n1, n2;
		n1 = new double[2];
		n2 = new double[2];

		resolDeg = (int) ( Math.round( ( 2*Math.asin(1f/(2 * scale))*(1f/4) / TwoPI) * 360 ) );
		resolDeg = (resolDeg>=1)? resolDeg : 1;

		double r1 = neuronDiam*scale-neuronDiam/2;
		double r2 = neuronDiam*scale+neuronDiam/2;// neuronDiam length

		offsets = new ArrayList<ArrayList<double[]>>();

		for (int aDeg = 0; aDeg<360; aDeg+=resolDeg) {

			float aRad = ((float)aDeg/360)*TwoPI;

			ArrayList<double[]> offsetsAngle = new ArrayList<double[]>();

			n1[0] = r1 * (float) Math.cos(aRad);
			n1[1] = r1 * (-1) * (float) Math.sin(aRad);

			n2[0] = r2 * (float) Math.cos(aRad);
			n2[1] = r2 * (-1) * (float) Math.sin(aRad);

//			IJ.log(""+ n1[0] + " , " + n1[1]);
//			IJ.log(""+ n2[0] + " , " + n2[1]);

			for (double x = -outerRange; x <= outerRange; x+=samplingStep) {
				for (double y = -outerRange; y <= outerRange; y+=samplingStep) {

					double px =  x;//x-size/2+0.5;
					double py = -y;//-(y-size/2+0.5);

					double dst = point2line(n1[0], n1[1], n2[0], n2[1], px, py);

					if (dst<=neuronDiam/2) {
						offsetsAngle.add(new double[]{px, py});
					}

				}
			}

			offsets.add(offsetsAngle);

		}

		/*
		 set profiles
		*/
		profiles = new float [locations.length][offsets.size()];

//        for (int i=0; i<1; i++) {
//			IJ.log("start"+offsets.get(i).size());
//            for (int j = 0; j<offsets.get(i).size(); j++) {
//                IJ.log(i+", "+offsets.get(i).get(j)[0]+", "+offsets.get(i).get(j)[1]+"\n");
//            }
//        }

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
                                );
                    }

                    profiles[locIdxProfile][offsetIdx] /= offsets.get(offsetIdx).size();

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

		double proj = (p_b[0] - n2x) * (n1x-n2x) + (p_b[1] - n2y) * (n1y-n2y);

		if(proj<0){
			return Double.MAX_VALUE;
		}

		proj = (p_b[0] - n1x) * n[0] + (p_b[1] - n1y) * n[1];
		if(proj<0){
			return Double.MAX_VALUE;
		}

		//IJ.log("nLen: "+nLen+" -> "+n[0]+","+n[1]+" proj: "+proj);

		p_b[0] = p_b[0] - n1x - proj * n[0];
		p_b[1] = p_b[1] - n1y - proj * n[1];

//		distance_from_line = vectorNorm(distance_2d);

		return Math.sqrt(p_b[0]*p_b[0]+p_b[1]*p_b[1]);
	}


}
