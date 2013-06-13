package advantra.critpoint;

import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.filter.Convolver;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 5/25/13
 * Time: 3:37 PM
 * To change this template use File | Settings | File Templates.
 */

public class CircularConfiguration3 {

	public int			r;
	public int			d;

	public static int			nPeaks = 3;
	public static int 			nRot = 8;
	public static float 		rotStp = (float) ((2*Math.PI)/nRot);

	public static float TwoPi = (float) (2*Math.PI);

	ArrayList<float[][]> kernels = new ArrayList<float[][]>();
	ArrayList<String> names = new ArrayList<String>();

	public CircularConfiguration3(
		int 	patchRadius
	)
	{

		r = patchRadius;
		d = 2 * patchRadius + 1;

		Vector<int[]> comb = new Vector<int[]>();

		for (int Rpix = r; Rpix >= (int)(1.0*r); Rpix*= 0.5) {

			//System.out.println("start from "+(int) Math.round(0.6*Rpix));
			for (float ratio = 0.5f; ratio >= 0.5; ratio*=0.75){

				int Tpix = (int) Math.round(ratio*Rpix);

                if (Tpix>=2){
					int angRes = (int) Math.round( (ratio/TwoPi)*360 );
					angRes = (angRes/10 + 1)*10;
					comb.add(new int[]{Rpix, Tpix, angRes}); System.out.println("3 : "+Rpix+","+Tpix+","+angRes);
				}
				else {
					break;
				}
			}

		}

        Vector<int[]> comb1 = new Vector<int[]>();

		for (int combIdx = 0; combIdx < comb.size(); combIdx++) {

			int minAngScaleDegrees = comb.get(combIdx)[2];
			int Rpx = comb.get(combIdx)[0];
			int Tpx = comb.get(combIdx)[1];

			int startscale = 2;

			for (int d1 = startscale*minAngScaleDegrees; d1 < 360; d1+=minAngScaleDegrees){
				for (int d2 = startscale*minAngScaleDegrees; d2 < 360; d2+=minAngScaleDegrees){
					for (int d3 = startscale*minAngScaleDegrees; d3 < 360; d3+=minAngScaleDegrees){

						boolean isConfiguration = (d1+d2+d3==360) && (d1<=180) && (d2<=180) && (d3<=180);

						if(isConfiguration) {

							boolean covered = false;
							// check if it exists so far in other rotations
							for (int k = 0; k < comb1.size(); k++){
								if(
									(
													d1        ==comb1.get(k)[0] &&
													d2        ==comb1.get(k)[1] &&
													d3        ==comb1.get(k)[2] &&
													Rpx		  ==comb1.get(k)[3] &&
													Tpx		  ==comb1.get(k)[4]
									)
									||
									(
													d1        ==comb1.get(k)[1] &&
													d2        ==comb1.get(k)[2] &&
													d3        ==comb1.get(k)[0] &&
													Rpx		  ==comb1.get(k)[3] &&
													Tpx       ==comb1.get(k)[4]
									)
									||
									(
													d1        ==comb1.get(k)[2] &&
													d2        ==comb1.get(k)[0] &&
													d3        ==comb1.get(k)[1] &&
													Rpx		  ==comb1.get(k)[3] &&
													Tpx       ==comb1.get(k)[4]
									)
								)
								{
									covered = true;
								}
							}

							if(!covered){

								comb1.add(new int[]{d1, d2, d3, Rpx, Tpx});
								float[] inhere = new float[nPeaks];
								inhere[0] = (d1/360f)*TwoPi;
								inhere[1] = (d2/360f)*TwoPi;
								inhere[2] = (d3/360f)*TwoPi;
								String name = ""+d1+","+d2+","+d3+","+Rpx+","+Tpx;
								names.add(name);
								//System.out.println("3 forming : "+name);
								kernels.add(formKernel(inhere, Rpx, Tpx));

							}

						}

					}

				}
			}
		}
	}

	public ImageStack plotKernels()
	{

		ImageStack is = new ImageStack(d, d);

		for (int i = 0; i < kernels.size(); i++) {
			 ImageProcessor ip = new FloatProcessor(d, d, kernels.get(i)[0]);
			is.addSlice(names.get(i), ip);
		}

		return is;

	}

	public ImageStack plotKernel(
		int kerIdx
	)
	{

		ImageStack is = new ImageStack(d, d);

		for (int i = 0; i < nRot; i++) {
			ImageProcessor ip = new FloatProcessor(d, d, kernels.get(kerIdx)[i]);
			is.addSlice(names.get(kerIdx), ip);
		}

		return is;

	}

	public ImageProcessor plotKernel(
		int kerIdx,
		int rotIdx
	)
	{

		ImageProcessor ip = new FloatProcessor(d, d, kernels.get(kerIdx)[rotIdx]);

		return ip;

	}

	/*
	SCORE CALCULATION WHOLE IMAGE
	 */

	public ImageProcessor score(
		int kernelIdx,
		ImageProcessor input
	)
	{
		ImageProcessor ip = score(kernelIdx, 0, input);  // rot 0

		for (int i = 1; i < nRot; i++) {
			ImageProcessor ip_sc = score(kernelIdx, i, input);
			for (int j = 0; j < ip_sc.getHeight()*ip_sc.getWidth(); j++) {
				if (ip_sc.getf(j) > ip.getf(j)) {
					ip.setf(j, ip_sc.getf(j));
				}
			}
		}

		return ip;

	}

	public ImageStack scoreAllRot(
		int kernelIdx,
		ImageProcessor input
	)
	{
	  	ImageStack is = new ImageStack(input.getWidth(), input.getHeight());

		for (int i = 0; i < nRot; i++) {
			is.addSlice(score(kernelIdx, i, input));
		}

		return is;

	}

	public ImageProcessor score(
		int kernelIdx,
		int rotIdx,
		ImageProcessor input
	)
	{
		ImageProcessor ip = input.duplicate();//new FloatProcessor(input.getWidth(), input.getHeight(), (float[]) input.getPixels());

		//convolution
		Convolver c = new Convolver();
		c.setNormalize(false); // important not to normalize (did my own normalization)

		float[] k = kernels.get(kernelIdx)[rotIdx];
		//new ImagePlus("before."+rotIdx, input.duplicate()).show();
		//new ImagePlus("kernel."+rotIdx, new FloatProcessor(d, d, k)).show();
		c.convolveFloat(ip, k, d, d);
		//new ImagePlus("after."+rotIdx, ip).show();

		return ip;
	}

    public ImageStack score_Experimental(
            int kernelIdx,
            int rotIdx,
            ImageProcessor input
//            int atX,
//            int atY
    )
    {

//        ImageProcessor ip = input.duplicate();//new FloatProcessor(input.getWidth(), input.getHeight(), (float[]) input.getPixels());
//        //convolution
//        Convolver c = new Convolver();
//        c.setNormalize(false); // important not to normalize (did my own normalization)

        ImageProcessor ipOut = new FloatProcessor(input.getWidth(), input.getHeight());

        float[] k = kernels.get(kernelIdx)[rotIdx];

        float val = 0;

        for (int imgX = r; imgX < input.getWidth()-r; imgX++) {
            for (int imgY = r; imgY < input.getHeight()-r; imgY++) {

                System.out.print("x: "+imgX+" , y: "+imgY+"    ");

                int cntP = 0;
                int cntN = 0;

                // Welford algorithm

                for (int ki = 0; ki < k.length; ki++) {  // take mean and variance estimate in one loop

                    if (k[ki]>0) {
                        cntP++;
                    }
                    else if(k[ki]<0) {
                        cntN++;
                    }


                }

                System.out.println("positives counted:  "+cntP+" , negatives counted "+cntN);
            }
        }

        return null;


        //new ImagePlus("before."+rotIdx, input.duplicate()).show();
        //new ImagePlus("kernel."+rotIdx, new FloatProcessor(d, d, k)).show();
        //c.convolveFloat(ip, k, d, d);
        //new ImagePlus("after."+rotIdx, ip).show();
        //return ip;
    }

	/*
	SCORE CALCULATION POSITION
	 */

	public float score(
		int atX,
		int atY,
		int kernelIdx,
		ImageProcessor input
	)
	{
		float scMax = Float.NEGATIVE_INFINITY;

		for (int rot = 0; rot < nRot; rot++) {

			float sc = 0;

			for (int x = 0; x < d; x++) {
				for (int y = 0; y < d; y++) {

					int x_img = x + atX - r;
					int y_img = y + atY - r;

					if (x_img>=0 && x_img<input.getWidth() && y_img>=0 && y_img<input.getHeight()) {
						sc += kernels.get(kernelIdx)[rot][x+d*y] * input.getf(x_img, y_img);
					}

				}
			}

			if(sc>scMax) {
				scMax = sc;
			}

		}

		return scMax;

	}

    /*
    AUX METHODS
     */

	private float[][] formKernel(
		float[] angles,
		int Rpix,
		int Tpix
	)
	{

		float[][] peaksRad = new float[nRot][nPeaks];

		for (int cnt_rots = 0; cnt_rots < nRot; cnt_rots++) {

			float start_pos = cnt_rots*rotStp;

			for (int cnt_pks = 0; cnt_pks < nPeaks; cnt_pks++) {

				if(cnt_pks==0)
					peaksRad[cnt_rots][cnt_pks] 	= start_pos;
				else
					peaksRad[cnt_rots][cnt_pks] 	=
							peaksRad[cnt_rots][cnt_pks-1] +
									angles[cnt_pks-1];

			}

		}

		float[][] kernels = new float[nRot][];

		int kernelsW = d;
		int kernelsH = d;

		// form the kernels
		int xc = kernelsW/2;
		int yc = kernelsH/2;
		int[] cent = new int[]{xc, yc};
		double[] n = new double[2];
		int[] p = new int[2];

		for (int rIdx = 0; rIdx < nRot; rIdx++) {

			kernels[rIdx] = new float[kernelsH*kernelsW];

			int nrON = 0, nrOFF = 0;

			for (int x = 0; x < kernelsW; x++) {
				for (int y = 0; y < kernelsH; y++) {
					if ( (x-xc)*(x-xc)+(y-yc)*(y-yc) <= Rpix*Rpix){

						boolean isON = false;

						for (int pIdx = 0; pIdx < nPeaks; pIdx++) {
							double ang = peaksRad[rIdx][pIdx];

							n[0] = Math.sin(ang);
							n[1] = -Math.cos(ang);

							p[0] = x;
							p[1] = y;

							double dst = distance_point_to_line_2d(cent, n, p);
							if (Math.round(dst) < Tpix/2) { // <= Tpix/2
								kernels[rIdx][x+kernelsW*y] = 1;
								nrON++;
								isON = true;
								break;
							}

						}

						if(!isON) {

							kernels[rIdx][x+kernelsW*y] = -1;
							nrOFF++;

						}

					}
					else {
						kernels[rIdx][x+kernelsW*y] = 0;
					}
				}
			}

			// normalize
			for (int x = 0; x < kernelsW; x++) {
				for (int y = 0; y < kernelsH; y++) {
					if (kernels[rIdx][x+kernelsW*y] > 0) {
						kernels[rIdx][x+kernelsW*y] /= nrON;
					}
					else if (kernels[rIdx][x+kernelsW*y] < 0) {
						kernels[rIdx][x+kernelsW*y] /= nrOFF;
					}
				}
			}

		}

		return kernels;

	}

	private double 	distance_point_to_line_2d(
            int[] line_base_point,
            double[] line_unit_direction,
            int[] point
    )
    {
		// line is in vector form (point+unit_direction)

		double distance_from_line = 0;

		double[] distance_2d = new double[2];

		// line_base_point - point
		distance_2d[0] = line_base_point[0] - point[0];
		distance_2d[1] = line_base_point[1] - point[1];
//        distance_2d[2] = line_base_point[2] - point[2]; 	// it contains difference a-p

		double proj =
				distance_2d[0] * line_unit_direction[0] +
				distance_2d[1] * line_unit_direction[1];// +
//                        distance_2d[2] * line_unit_direction[2]; 	// dotproduct((a-p), n)

		if(proj<0){
			return vectorNorm(distance_2d);
		}

		distance_2d[0] = distance_2d[0] - proj * line_unit_direction[0];
		distance_2d[1] = distance_2d[1] - proj * line_unit_direction[1];
//        distance_2d[2] = distance_2d[2] - proj * line_unit_direction[2]; // || (a-p) - ((a-p)*n) * n ||

		distance_from_line = vectorNorm(distance_2d);

		return distance_from_line;
	}

	private double 	vectorNorm(
		double[] vector
	)
	{
		double norm = 0;

		for (int i = 0; i < vector.length; i++) {
			norm += vector[i] * vector[i];
		}

		return Math.sqrt(norm);
	}

}