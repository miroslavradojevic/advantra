package advantra.critpoint;

import ij.IJ;
import ij.ImageStack;
import ij.plugin.filter.Convolver;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 5/29/13
 * Time: 3:37 AM
 * To change this template use File | Settings | File Templates.
 */
public class SymmConfiguration {

	public int			r;
	public int			d;

	public static int			maxNrOnRegions = 1;
	public static int 			nRot = 15;
	public static float 		rotStp = (float) ((2*Math.PI)/nRot);

	public static float TwoPi = (float) (2*Math.PI);

	ArrayList<float[][]> kernels = new ArrayList<float[][]>();
	ArrayList<String> names = new ArrayList<String>();

	public SymmConfiguration(
		int 	patchRadius
	)
	{
		r = patchRadius;
		d = 2 * patchRadius + 1;

		Vector<int[]> comb = new Vector<int[]>();

		for (int Rpix = r; Rpix >= (int)(0.5*r); Rpix*= 0.75) {

			for (int nrOnRegions = 1; nrOnRegions <= maxNrOnRegions; nrOnRegions+=1){
				comb.add(new int[]{Rpix, nrOnRegions});
			}

		}

		for (int combIdx = 0; combIdx < comb.size(); combIdx++) {

			int Rpx = comb.get(combIdx)[0];
			int Nrg = comb.get(combIdx)[1];

			float[][] inhere = new float[Nrg][2];

			String name = Rpx+",";

			for (int i = 0; i < Nrg; i++) {
				inhere[i][0] = (2*i)*(TwoPi/(2*Nrg)); inhere[i][1] = (2*i+1)*(TwoPi/(2*Nrg));
				name+="["+IJ.d2s(inhere[i][0],2)+","+IJ.d2s(inhere[i][1],2)+"]";
			}

			names.add(name);
			kernels.add(formKernel(inhere, Rpx));

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

		c.convolveFloat(ip, k, d, d);

		return ip;

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
										float[][] intervals,
										int Rpix
	)
	{

		float[][][] intervalsRot = new float[nRot][intervals.length][];

		for (int cnt_rots = 0; cnt_rots < nRot; cnt_rots++) {

			float start_pos = 0;//*rotStp;

			for (int cnt_itvs = 0; cnt_itvs < intervals.length; cnt_itvs++) {

				intervalsRot[cnt_rots][cnt_itvs] 	= new float[]{cnt_rots*rotStp+intervals[cnt_itvs][0], cnt_rots*rotStp+intervals[cnt_itvs][1]};

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

						for (int intvIdx = 0; intvIdx < intervals.length; intvIdx++) {

							float ang1 = wrap_0_2PI(intervalsRot[rIdx][intvIdx][0]);
							float ang2 = wrap_0_2PI(intervalsRot[rIdx][intvIdx][1]);

							float ang = wrap_0_2PI((float) (Math.atan2(y-yc, x-xc)+Math.PI));

//							n[0] = Math.sin(ang);
//							n[1] = -Math.cos(ang);
//
//							p[0] = x;
//							p[1] = y;

//							double dst = distance_point_to_line_2d(cent, n, p);
							if (wrap_PI(ang - ang1) >= 0 && wrap_PI(ang-ang2) <= 0) {
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


	private float wrap_0_2PI(
		float in
	)
	{

		float out = in;

		while(out<0){
			out += 2*Math.PI;
		}
		while(out>=2*Math.PI){
			out -= 2*Math.PI;
		}

		return out;
	}

	private static double wrap_PI(
		float in
	)
	{

		float out = in;

		while(out<=-Math.PI){
			out += 2*Math.PI;
		}
		while(out>Math.PI){
			out -= 2*Math.PI;
		}

		return out;
	}

}
