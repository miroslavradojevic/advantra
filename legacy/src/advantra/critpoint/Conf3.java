package advantra.critpoint;

import ij.IJ;
import ij.ImageStack;
import ij.gui.Plot;
import ij.plugin.filter.Convolver;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 5/31/13
 * Time: 2:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class Conf3 {

	public int			r;
	public int			d;

	public static int			nPeaks = 3;
	public static int 			nRot = 12;
	public static float 		rotStp = (float) ((2*Math.PI)/nRot);

	public static float TwoPi = (float) (2*Math.PI);

	float sCenter;
	int nCenter;

	ArrayList<int[][]> regionIdxMap = new ArrayList<int[][]>(); // each element is a template with indexes for each rotation
	ArrayList<int[][]>   regionSize	= new ArrayList<int[][]>();
	ArrayList<String> names = new ArrayList<String>();

	// aim is to calculate average in each region (2*nrPeaks+1 regions)
	float[] regionSums; // this variable will store sums for current rotation regions when calculating score at some location
	float Xon, Xoff;  // parameters to set the score for each region

	public Conf3(
		int 	patchRadius
	)
	{

		r = patchRadius;
		r = (r<3)? 3 : r; // lower limit
		d = 2*r+1;

		regionSums = new float[2*nPeaks+1];
		Xon = 200;
		Xoff = 50;

		//create different combinations radius-thickness
		Vector<int[]> comb = new Vector<int[]>();

		for (int R2 = r; R2 > (int)(0.5*r); R2*= 0.75) {

			for (float ratio = 0.5f; ratio > 0.1; ratio*=0.75){

				int T = (int) Math.round(ratio*R2);

				if (T>=2){

					int angRes = (int) Math.round( (ratio/TwoPi)*360 );
					angRes = (angRes/20 + 1)*20;

					// add configuration
					/*
						R2 - ring radius
					    R1 - ring radius
					    R3  - inner ring radius
					*/
					int R1 = (int)Math.round(0.5*R2);
					R1 = (R1<=2)? 3 : R1;
					int R3 = (int)Math.round(0.2*R2);
					R3 = (R3==R1)? R1-1 : R3;

					comb.add(new int[]{R2, R1, R3, T, angRes});

				}
				else {
					break;
				}

			}

		}

		// create different angular combinations
		Vector<int[]> comb1 = new Vector<int[]>();

		for (int combIdx = 0; combIdx < comb.size(); combIdx++) {

			int Rpx = comb.get(combIdx)[0];
			int Rlower = comb.get(combIdx)[1];
			int Rinner = comb.get(combIdx)[2];
			int Tpx = comb.get(combIdx)[3];
			int minAngScaleDegrees = comb.get(combIdx)[4];

			int startscale = 1;

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

								int[][] idxMap = new int[nRot][]; // nRot x (d*d)
								int[][] regSzs = new int[nRot][]; // nRot x (2*nPeaks+1)

								boolean isOK = formIdxMap(inhere, Rpx, Rlower, Rinner, Tpx, idxMap, regSzs);

								if(isOK) {

									String name = ""+d1+","+d2+","+d3+","+Rpx+","+Tpx;
									names.add(name);
									regionIdxMap.add(idxMap);
									regionSize.add(regSzs);

								}

							}

						}

					}

				}
			}
		}

		System.out.println("formed "+regionIdxMap.size()+" configurations.");

	}

	public ImageStack plotAll()
	{

		ImageStack is = new ImageStack(d, d);

		for (int i = 0; i < regionIdxMap.size(); i++) {
			ImageProcessor ip = new FloatProcessor(d, d, regionIdxMap.get(i)[0]);
			is.addSlice(names.get(i), ip);
		}

		return is;

	}

	public ImageStack plotOne(
		int kernelIdx
	)
	{

		ImageStack is = new ImageStack(d, d);

		for (int i = 0; i < nRot; i++) {
            System.out.println("rot "+i+" -> "+regionIdxMap.get(kernelIdx)[i].length);
			ImageProcessor ip = new FloatProcessor(d, d, regionIdxMap.get(kernelIdx)[i]);
			is.addSlice("rot."+i+","+names.get(i), ip);
		}

		return is;

	}

	/*
	SCORE CALCULATION WHOLE IMAGE // TODO: finish
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

	public ImageProcessor score(// TODO: here first!
									   int kernelIdx,
									   int rotIdx,
									   ImageProcessor input
	)
	{
		ImageProcessor ip = input.duplicate();//new FloatProcessor(input.getWidth(), input.getHeight(), (float[]) input.getPixels());

		//convolution
		Convolver c = new Convolver();
		c.setNormalize(false); // important not to normalize (did my own normalization)

		//float[] k = kernels.get(kernelIdx)[rotIdx];

		//new ImagePlus("before."+rotIdx, input.duplicate()).show();
		//new ImagePlus("kernel."+rotIdx, new FloatProcessor(d, d, k)).show();

		//c.convolveFloat(ip, k, d, d);

		//new ImagePlus("after."+rotIdx, ip).show();

		return ip;

	}

	/*
	SCORE CALCULATION POSITION
	 */

	public double scoreAtPos(
		int atX,
		int atY,
		int kernelIdx,
		ImageProcessor input,
        float Xon,
        float Xoff
	)
	{
		double highestScr = 0;
        int getIndex;   // variable to store the index each time

		for (int rot = 0; rot < nRot; rot++) {

			//reset region sums for each rotation each time (max. will be taken)
            for (int t = 0; t < regionSums.length; t++) {
                regionSums[t] = 0;
            }

            for (int l = 0; l < d*d; l++) { // check all template locations for sum indexes
                getIndex  =  regionIdxMap.get(kernelIdx)[rot][l];
                if (getIndex>=0 && getIndex<=2*nPeaks) { // if the index is meaningful

                    int x = l % d; // col
                    int y = l / d; // row

                    int x_img = x + atX - r;
                    int y_img = y + atY - r;

                    if (x_img>=0 && x_img<input.getWidth() && y_img>=0 && y_img<input.getHeight()) { // if value can be read from this location
                        regionSums[ regionIdxMap.get(kernelIdx)[rot][l] ]+=input.getf(x_img, y_img);
                    }
                    else {
                        regionSums[ regionIdxMap.get(kernelIdx)[rot][l] ]+=0;
                    }
                }
            }

            double score = 1;

            //System.out.println(" before score... rotation "+rot);

            float[] x_axis_test = new float[regionSums.length];
            float[] lhoods = new float[regionSums.length];
            for (int t = 0; t < regionSums.length; t++) {
                x_axis_test[t] = t;
            }

            for (int t = 0; t < regionSums.length; t++) {

                //System.out.println("region "+t+" sum =  "+regionSums[t]+"  , size = "+regionSize.get(kernelIdx)[rot][t]+" avg: ");

                regionSums[t] /= regionSize.get(kernelIdx)[rot][t]; // t corresponds to region size at this rotation

                //System.out.print(t + " -> " + regionSums[t] + "  ,   ");

                // now make the final calculation
                if (t==0) { // it is central one
                    score *= 1 - Math.exp(-regionSums[0]/Xon); // higher value -> 1
                    //System.out.println(1 - Math.exp(-regionSums[0]/Xon));
                    lhoods[t] = (float) (1 - Math.exp(-regionSums[t]/Xon));
                }
                if (t>=1 && t <= nPeaks) { // peak regions
                    score *= 1 - Math.exp(-regionSums[t]/Xon);
                    //System.out.println(1 - Math.exp(-regionSums[t]/Xon));
                    lhoods[t] = (float) (1 - Math.exp(-regionSums[t]/Xon));
                }
                if (t>=nPeaks+1 && t <= 2*nPeaks) {
                    score *= Math.exp(-regionSums[t]/Xoff);
                    //System.out.println(Math.exp(-regionSums[t]/Xoff));
                    lhoods[t] = (float) (Math.exp(-regionSums[t]/Xoff));
                }
            }

            /*Plot p = new Plot("rot"+rot, "region", "avg", x_axis_test, regionSums);
            p.show();

            p = new Plot("rot"+rot, "region", "lhood", x_axis_test, lhoods);
            p.show();
            */

			if(score>highestScr) {
                highestScr = score;
			}

		}

		return highestScr;

	}

	/*
    AUX METHODS
     */

	private boolean formIdxMap(
		float[] 	angles,
		int 		Rpix,
		int 		Rlower,
		int 		Rinner,
		int 		Tpix,
		int[][] 	idxMap,
		int[][] 	regSzs
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

		// form the kernels
		int xc = d/2;
		int yc = d/2;
		int[] cent = new int[]{0, 0};
		float[] n = new float[2];
		int[] p = new int[2];
		float[] ap = new float[nPeaks];

		for (int rIdx = 0; rIdx < nRot; rIdx++) {

			idxMap[rIdx] = new int[d*d];
			regSzs[rIdx] = new int[2*nPeaks+1];

			for (int pIdx = 0; pIdx < nPeaks; pIdx++)
				ap[pIdx] = wrap_0_2PI( peaksRad[rIdx][pIdx] );

			Arrays.sort(ap); // wrapped and sorted angles for this rotation

//			int nrON = 0, nrOFF = 0;   // count those that are expected positive

			for (int x = 0; x < d; x++) {
				for (int y = 0; y < d; y++) {

					int d2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);

					if (d2 <= Rpix*Rpix && d2 >= Rlower*Rlower) {

						p[0] = x-xc;
						p[1] = -(y-yc);

						boolean isON = false;

						for (int pIdx = 0; pIdx < nPeaks; pIdx++) {

							float ang = peaksRad[rIdx][pIdx];  		// read the angle

							n[0] = (float) Math.cos(ang);
							n[1] = (float) Math.sin(ang);

							float dst = point2dir(cent, n, p);

							if ( Math.round(dst) <= Tpix/2) { // belongs to pIdx ON peak

								idxMap[rIdx][x+d*y] = pIdx+1;
								regSzs[rIdx][pIdx+1]++;
								isON = true;
								break;

							}

						}

						if(!isON) {

							float a = wrap_0_2PI((float) Math.atan2(p[1], p[0]));

							boolean added = false;
							for (int pIdx = 0; pIdx < nPeaks-1; pIdx++) {

								if (a>ap[pIdx] && a<=ap[pIdx+1]) {
									idxMap[rIdx][x+d*y] = nPeaks+1+pIdx;
									regSzs[rIdx][nPeaks+1+pIdx]++;
									added = true;
									break;
								}

							}
							if (!added) {
								idxMap[rIdx][x+d*y] = 2*nPeaks;
								regSzs[rIdx][2*nPeaks]++;
							}

						}

					}
					else if ( d2 <= Rinner*Rinner ) {
						idxMap[rIdx][x+d*y] = 0;
						regSzs[rIdx][0]++;
					}
					else {
						idxMap[rIdx][x+d*y] = -1; // nothing
					}
				}
			}

			//System.out.println("finished defining the template for rotation idx. "+rIdx+" : ");

			for (int k = 0; k < 2*nPeaks+1; k++) {

				if (regSzs[rIdx][k] <= 3 ) {
					//System.out.print("\n configuration: \n");
					//System.out.print("R2="+Rpix+",R1="+Rlower+",R3="+Rinner+",T="+Tpix+"\n");
					for (int l = 0; l < angles.length; l++) {
						System.out.print("a("+l+") = "+((angles[l]/TwoPi)*360));
					}
					System.out.print("\n <= 3 samples in region "+k+" at rot. "+rIdx);
					return false;
				}

			}

		}

		return true;

	}

	private float 	point2dir(
		int[] 	b,    	// direciton base point
		float[] n, 		// direction vector
		int[] 	p     	// point considered
	)
	{
		// line is in vector from b in n direction

		float d = 0;

		float[] p_b = new float[2];

		// p - b
		p_b[0] = p[0] - b[0];
		p_b[1] = p[1] - b[1];

		float proj = p_b[0] * n[0] + p_b[1] * n[1];

		if(proj<0){
			// "behind" the orientation
			return (float) Math.sqrt(p_b[0]*p_b[0]+p_b[1]*p_b[1]);
		}

		// || (p-b) - dot(p-b,n) * n ||
		p_b[0] = p_b[0] - proj * n[0];
		p_b[1] = p_b[1] - proj * n[1];

//		distance_from_line = vectorNorm(distance_2d);

		return (float) Math.sqrt(p_b[0]*p_b[0]+p_b[1]*p_b[1]);
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

	private static float wrap_PI(
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
