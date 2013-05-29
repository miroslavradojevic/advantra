package advantra.critpoint;

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

	public static int			maxNrOnRegions = 2;
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

//				int Tpix = (int) Math.round(ratio*Rpix);
//
//				if (Tpix>=2){
//					int angRes = (int) Math.round( (ratio/TwoPi)*360 );
//					angRes = (angRes/20 + 1)*20;
//					comb.add(new int[]{Rpix, Tpix, angRes});
					comb.add(new int[]{Rpix, nrOnRegions});
			//kernels.add(formKernel(new float[]{0, TwoPi/nBorders}, Rpix));
//				}
//				else {
//					break;
//				}
			}

		}

		//System.out.println(""+comb.size()+" R combinations formed");

		for (int combIdx = 0; combIdx < comb.size(); combIdx++) {

			int Rpx = comb.get(combIdx)[0];
			int Nrg = comb.get(combIdx)[1];

			float[][] inhere = new float[Nrg][2];

			for (int i = 0; i < Nrg; i++) {
				inhere[i][0] = (2*i)*(TwoPi/(2*Nrg)); inhere[i][1] = (2*i+1)*(TwoPi/(2*Nrg));
			}

			kernels.add(formKernel(inhere, Rpx));

		}
	}

	private float[][] formKernel(
										float[][] intervals,
										int Rpix
	)
	{

		float[][][] intervalsRot = new float[nRot][intervals.length][];

		for (int cnt_rots = 0; cnt_rots < nRot; cnt_rots++) {

			float start_pos = 0;//*rotStp;

			for (int cnt_itvs = 0; cnt_itvs < intervals.length; cnt_itvs++) {

//				if(cnt_itvs==0)
				intervalsRot[cnt_rots][cnt_itvs] 	= new float[]{cnt_rots*nRot+intervals[cnt_itvs][0], cnt_rots*nRot+intervals[cnt_itvs][1]};
//				else
//					intervalsRot[cnt_rots][cnt_itvs] 	= new float[]{start_pos+cnt_rots*rotStp+intervals[cnt_itvs][0], start_pos+cnt_rots*rotStp+intervals[cnt_itvs][1]};
//							peaksRad[cnt_rots][cnt_pks-1] +
//									angles[cnt_pks-1];

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
							double ang1 = intervalsRot[rIdx][intvIdx][0];
							double ang2 = intervalsRot[rIdx][intvIdx][1];

							double ang = Math.atan2(y-yc, x-xc)+Math.PI;

//							n[0] = Math.sin(ang);
//							n[1] = -Math.cos(ang);
//
//							p[0] = x;
//							p[1] = y;

//							double dst = distance_point_to_line_2d(cent, n, p);
							if (ang >= ang1 && ang <= ang2) {
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

}
