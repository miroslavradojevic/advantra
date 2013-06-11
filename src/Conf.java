import ij.IJ;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/11/13
 * Time: 1:44 PM
 */
public class Conf {

    public int			r;
    public int			d;
    public int          diam;

    public static int   nPeaks = 3;
    public static int 	nRot = 8;
    public static float rotStp = (float) ((2*Math.PI)/nRot);
    public static float TwoPi = (float) (2*Math.PI);

    ArrayList<int[]> angles         = new ArrayList<int[]>();
    ArrayList<int[]> regionIdxMap   = new ArrayList<int[]>();
    ArrayList<int[]> regionSize	    = new ArrayList<int[]>();
    ArrayList<String> names         = new ArrayList<String>();
    ArrayList<int[]> kernels      = new ArrayList<int[]>(); // for each index map (configuration) -> set of rotations

    public Conf(int confRadius, int diam) {

        r = confRadius;
        r = (r<5)? 5 : r; // lower limit
        d = 2*r+1;

        this.diam = diam;

        // diam will define min angle step
        double  angStep     = 2 * Math.asin((float)diam/(2*r));
        int resol = 10;
        int     angStepDeg  = (int) Math.round( ((angStep/(2*Math.PI))*360) / resol) * resol;
        angStepDeg = (angStepDeg<resol)?resol:angStepDeg;
        System.out.println(angStep+" , step deg.  -> "+angStepDeg+" "+(((angStep/(2*Math.PI))*360)));
        //create different angular combinations of 3 angles
        // making 360 together,
        for (int a1 = angStepDeg; a1<360; a1+=resol) {

            for (int a2 = angStepDeg; a2<360; a2+=resol) {

                for (int a3 = angStepDeg; a3<360; a3+=resol) {

                    if(a1+a2+a3==360 && a1<=180 && a2<=180 && a3<=180) {

                        float[] angRad = new float[nPeaks];
                        angRad[0] = (a1/360f)*TwoPi;
                        angRad[1] = (a2/360f)*TwoPi;
                        angRad[2] = (a3/360f)*TwoPi;

                        //int[] idxMap1 = new int[d*d];     // (d*d)
                        //int[] regSzs1 = new int[2*nPeaks+1];
                        int[][] idxSzs = new int[2][];
                        //idxSzs[0] = idxMap1;
                       // idxSzs[1] = regSzs1;

                        boolean isOK = formIdxMap(angRad, r, 4, 3, diam, idxSzs);

                        if(isOK) {

                            // check whether it exists already
                            boolean covered = false;
                            // check if it exists so far in other rotations
                            for (int k = 0; k < angles.size(); k++){
                                if(
                                        (
                                                        a1        ==angles.get(k)[0] &&
                                                        a2        ==angles.get(k)[1] &&
                                                        a3        ==angles.get(k)[2]
                                        )
                                                ||
                                                (
                                                        a1        ==angles.get(k)[1] &&
                                                        a2        ==angles.get(k)[2] &&
                                                        a3        ==angles.get(k)[0]
                                                )
                                                ||
                                                (
                                                        a1        ==angles.get(k)[2] &&
                                                        a2        ==angles.get(k)[0] &&
                                                        a3        ==angles.get(k)[1]
                                                )
                                        )
                                {
                                    covered = true;
                                }
                            }

                            if (!covered) {
                                String name = ""+a1+","+a2+","+a3+","+r+","+diam;
                                names.add(name);
                                angles.add(new int[]{a1, a2, a3});
                                regionIdxMap.add(idxSzs[0]);
                                regionSize.add(idxSzs[1]);
                            }

                        }

                    }

                }

            }

        }

        System.out.println("total "+regionIdxMap.size()+" configurations!");

    }

    private boolean formIdxMap(
            float[] 	angles,
            int 		Rpix,
            int 		Rlower,
            int 		Rinner,
            float 		Tpix,
            int[][]   	AAA
    )
    {

        int[] idxMapLocal;
        int[] regSzsLocal;

        float[] peaksRad = new float[nPeaks];

        for (int cnt_pks = 0; cnt_pks < nPeaks; cnt_pks++) {

                if(cnt_pks==0)
                    peaksRad[cnt_pks] 	= 0;
                else
                    peaksRad[cnt_pks] 	=
                            peaksRad[cnt_pks-1] +
                                    angles[cnt_pks-1];

        }

//        for (int k = 0; k < nPeaks; k++)
//            System.out.println(k+" : "+peaksRad[k]);

        // form the kernels
        int xc = d/2;
        int yc = d/2;
        int[] cent = new int[]{0, 0};
        float[] n = new float[2];
        int[] p = new int[2];
        float[] ap = new float[nPeaks];

        idxMapLocal = new int[d*d];
        regSzsLocal = new int[2*nPeaks+1];

        for (int pIdx = 0; pIdx < nPeaks; pIdx++)
            ap[pIdx] = wrap_0_2PI( peaksRad[pIdx] );

        Arrays.sort(ap); // wrapped and sorted angles for this rotation

            for (int x = 0; x < d; x++) {
                for (int y = 0; y < d; y++) {

                    int d2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);

                    if (d2 <= Rpix*Rpix && d2 >= Rlower*Rlower) {

                        p[0] = x-xc;
                        p[1] = -(y-yc);

                        boolean isON = false;

                        for (int pIdx = 0; pIdx < nPeaks; pIdx++) {

                            float ang = peaksRad[pIdx];  		// read the angle

                            n[0] = (float) Math.cos(ang);
                            n[1] = (float) Math.sin(ang);

                            float dst = point2dir(cent, n, p);

                            if ( dst <= Tpix/2) { // belongs to pIdx ON peak

                                idxMapLocal[x+d*y] = pIdx+1;
                                regSzsLocal[pIdx+1]++;
                                isON = true;
                                break;

                            }

                        }

                        if(!isON) {

                            float a = wrap_0_2PI((float) Math.atan2(p[1], p[0]));

                            boolean added = false;
                            for (int pIdx = 0; pIdx < nPeaks-1; pIdx++) {

                                if (a>ap[pIdx] && a<=ap[pIdx+1]) {
                                    idxMapLocal[x+d*y] = nPeaks+1+pIdx;
                                    regSzsLocal[nPeaks+1+pIdx]++;
                                    added = true;
                                    break;
                                }

                            }
                            if (!added) {
                                idxMapLocal[x+d*y] = 2*nPeaks;
                                regSzsLocal[2*nPeaks]++;
                            }

                        }

                    }
                    else if ( d2 <= Rinner*Rinner ) {
                        idxMapLocal[x+d*y] = 0;
                        regSzsLocal[0]++;
                    }
                    else {
                        idxMapLocal[x+d*y] = -1; // nothing
                    }
                }
            }

        for (int k = 0; k < 2*nPeaks+1; k++) {

            if (regSzsLocal[k] <= 3 ) {
//                System.out.println("some reg was small "+k+"   : "+regSzsLocal[k]+" elements");
                return false;
            }

        }

        AAA[0] = idxMapLocal;
        AAA[1] = regSzsLocal;

//        System.out.println("AAA "+AAA.length+" x "+AAA[0].length);

        return true;

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

    /*
    plot
     */

    public ImageStack plot(
            int kernelIdx
    )
    {

        ImageStack is = new ImageStack(d, d);

        ImageProcessor ip = new FloatProcessor(d, d, regionIdxMap.get(kernelIdx));
        is.addSlice("kernel,"+names.get(kernelIdx), ip);

        return is;

    }


    public ImageStack plot()
    {

        System.out.println("index map : "+regionIdxMap.size());

        ImageStack is = new ImageStack(d, d);

        for (int l = 0; l<regionIdxMap.size(); l++) {
            ImageProcessor ip = new FloatProcessor(d, d, regionIdxMap.get(l));
            is.addSlice("kernel,"+names.get(l), ip);
        }

        return is;

    }

}