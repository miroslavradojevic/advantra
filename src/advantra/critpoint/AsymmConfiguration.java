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
 * Time: 8:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class AsymmConfiguration {

	public int			r;
	public int			d;

	public static int			N = 4;
	public static int 			nRot = 15;
	public static float 		rotStp = (float) ((2*Math.PI)/nRot);

	public static float TwoPi = (float) (2*Math.PI);

	ArrayList<float[][]> kernels = new ArrayList<float[][]>();
	ArrayList<String> names = new ArrayList<String>();

	public AsymmConfiguration(
		int 	patchRadius
	)
	{
		r = patchRadius;
		d = 2 * patchRadius + 1;

		Vector<int[]> comb = new Vector<int[]>();

		for (int Rpix = r; Rpix >= (int)(0.5*r); Rpix*= 0.75) {

			for (int i = 1; i < N; i++){

				int Tpix = (int)Math.round((i*Rpix)/(float)N);
				comb.add(new int[]{Rpix, Tpix});
			}

		}

		for (int combIdx = 0; combIdx < comb.size(); combIdx++) {

			int Rpx = comb.get(combIdx)[0];
			int Tpx = comb.get(combIdx)[1];

//			float[][] inhere = new float[Nrg][2];

			String name = Rpx+","+Tpx;

//			for (int i = 0; i < Nrg; i++) {
//				inhere[i][0] = (2*i)*(TwoPi/(2*Nrg)); inhere[i][1] = (2*i+1)*(TwoPi/(2*Nrg));
//				name+="["+ IJ.d2s(inhere[i][0], 2)+","+IJ.d2s(inhere[i][1],2)+"]";
//			}

			names.add(name);
			kernels.add(formKernel(Rpx, Tpx));

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
										int Rpix,
										int Tpix
	)
	{

		float[][] kernels = new float[nRot][];

		int kernelsW = d;
		int kernelsH = d;

		// form the kernels
		int xc = kernelsW/2;
		int yc = kernelsH/2;
		int[] cent = new int[]{xc, yc};
		double[] dir = new double[2];
		double[] p = new double[2];

		for (int rIdx = 0; rIdx < nRot; rIdx++) {  // nRot

			kernels[rIdx] = new float[kernelsH*kernelsW];

			int nrON = 0, nrOFF = 0;

			for (int x = 0; x < kernelsW; x++) {
				for (int y = 0; y < kernelsH; y++) {
					if ( (x-xc)*(x-xc)+(y-yc)*(y-yc) <= Rpix*Rpix ) {

//						boolean isON = false;

						float rotAngle = rIdx*rotStp;

						dir[0] = -Tpix * Math.sin(rotAngle);
						dir[1] = -Tpix * (-Math.cos(rotAngle));

						p[0] = x - (xc-dir[0]);
						p[1] = y - (yc-dir[1]);

						if ( p[0]*dir[0] + p[1]*dir[1] > 0 ) { //wrap_PI(ang - ang1) >= 0 && wrap_PI(ang-ang2) <= 0
							kernels[rIdx][x+kernelsW*y] = 1;
							nrON++;
//							isON = true;
//							break;
						}

//						if(!isON) {
						else {
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


}
