import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.filter.Convolver;
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

    public int          rInner;
    public int          rLower;

    public static int   nPeaks = 3;
    public static int 	nRot = 8;
    public static float rotStp = (float) ((2*Math.PI)/nRot);
    public static float TwoPi = (float) (2*Math.PI);

    ArrayList<int[]> angles         = new ArrayList<int[]>();  // as many as configurations
    ArrayList<String> names         = new ArrayList<String>();
    // with rotation
    ArrayList<int[]> regionIdxMap   = new ArrayList<int[]>();
    ArrayList<int[]> regionSize	    = new ArrayList<int[]>();
    ArrayList<float[]> kernels        = new ArrayList<float[]>(); // for each index map (configuration) -> set of rotations

    Convolver convolver;

    public Conf(int diam, double scale) {

        this.diam = diam;

        r = (int) (diam*scale);
        r = (r<6)? 6 : r; // lower limit
        d = 2*r+1;

        rInner = diam/2;
        rLower = r/2;

        convolver = new Convolver();
        convolver.setNormalize(false);

        // diam will define min angle step
        double  angStep     = 2 * Math.asin((float)diam/(2*r));
        int resol = 10;
        int     angStepDeg  = (int) Math.round( ((angStep/(2*Math.PI))*360) / resol) * resol;
        angStepDeg = (angStepDeg<resol)?resol:angStepDeg;

        //create different angular combinations of 3 angles
        // making 360 together,
        for (int a1 = angStepDeg; a1<360; a1+=resol) { // angStepDeg

            for (int a2 = angStepDeg; a2<360; a2+=resol) {   //

                for (int a3 = angStepDeg; a3<360; a3+=resol) {   //

                    if(a1+a2+a3==360 && a1<=180 && a2<=180 && a3<=180) {

                        float[] angRad = new float[nPeaks];
                        angRad[0] = (a1/360f)*TwoPi;
                        angRad[1] = (a2/360f)*TwoPi;
                        angRad[2] = (a3/360f)*TwoPi;

                        int[][][] index_map = new int[1][][];
                        int[][][] region_size = new int[1][][];
                        float[][][] kernel_rot = new float[1][][];

//                        System.out.println(""+a1+","+a2+","+a3+" : "+angStepDeg+" : "+resol);
//                        System.out.println(""+angRad[0]+","+angRad[1]+","+angRad[2]+" : "+angStepDeg+" : "+resol);

                        boolean isOK = createTemplates(angRad, r, rLower, rInner, diam, index_map, region_size, kernel_rot);

                        if (a1==0) {System.out.println("it was ok? "+isOK);}

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

                                for (int r = 0; r<nRot; r++) {
                                    regionIdxMap.add(index_map[0][r]);
                                    regionSize.add(region_size[0][r]);
                                    kernels.add(kernel_rot[0][r]);
                                }

                            }

                        }

                    }

                }

            }

        }

        System.out.println("added total "+regionIdxMap.size()+", "+angles.size()+" configurations, "+Conf.nRot+"  rotations each ");

    }

    private boolean createTemplates(
            float[] 	angles,
            int 		Rpix,
            int 		Rlower,
            int 		Rinner,
            float 		Tpix,
            // outputs
            int[][][]     index_map,
            int[][][]     region_size,
            float[][][]   kernels_rot
    )
    {

        int[][]     idxMapLocal = new int[nRot][];
        int[][]     regSzsLocal = new int[nRot][];
        float[][]   kernelLocal = new float[nRot][];

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

        System.out.println("forming with angles...");
        for (int cnt_rots = 0; cnt_rots < nRot; cnt_rots++) {
            for (int cnt_pks = 0; cnt_pks < nPeaks; cnt_pks++) {
                System.out.print(peaksRad[cnt_rots][cnt_pks]+", ");
            }
            System.out.println();
        }
        System.out.println("--\n");

        // form the kernels
        int xc = d/2;
        int yc = d/2;
        int[] cent = new int[]{0, 0};
        float[] n = new float[2];
        int[] p = new int[2];
        float[] ap = new float[nPeaks];

        for (int rIdx = 0; rIdx < nRot; rIdx++) {

            idxMapLocal[rIdx] = new int[d*d];
            regSzsLocal[rIdx] = new int[2*nPeaks+1];
            kernelLocal[rIdx] = new float[d*d];



            for (int pIdx = 0; pIdx < nPeaks; pIdx++)
                ap[pIdx] = ( peaksRad[rIdx][pIdx] );

            for (int i = 0; i<3; i++) {System.out.print("ap in:"+ap[i]);}

            Arrays.sort(ap);

            for (int i = 0; i<3; i++) {System.out.print("ap srt:"+ap[i]);}

            for (int pIdx = 0; pIdx < nPeaks; pIdx++)
                ap[pIdx] = wrap_0_2PI( peaksRad[rIdx][pIdx] );

             // wrapped and sorted angles for this rotation

            for (int i = 0; i<3; i++) {System.out.print("ap wrapped:"+ap[i]);}

            int sumON = 0;
            int sumOFF = 0;// for kernel normalization

            for (int x = 0; x < d; x++) {
                for (int y = 0; y < d; y++) {

                    int d2 = (x-xc)*(x-xc)+(y-yc)*(y-yc);

                    if (d2 <= Rpix*Rpix && d2 >= Rlower*Rlower) {

                        p[0] = x-xc;
                        p[1] = -(y-yc);

                        boolean isON = false;

                        for (int pIdx = 0; pIdx < nPeaks; pIdx++) {      // check first peak first always

                            float ang = peaksRad[rIdx][pIdx];  		// read the angle

                            n[0] = (float) Math.cos(ang);
                            n[1] = (float) Math.sin(ang);

                            float dst = point2dir(cent, n, p);

                            if (dst<=Tpix/2 && idxMapLocal[rIdx][x+d*y]==0) { // belongs to pIdx ON peak and not filled

                                idxMapLocal[rIdx][x+d*y] = pIdx+1;
                                regSzsLocal[rIdx][pIdx+1]++;
                                kernelLocal[rIdx][x+d*y] = +1;
                                sumON++;
                                isON = true;
                                //break; check all three peaks independently

                            }

                        }

                        if(!isON) {

                            float a = wrap_0_2PI((float) Math.atan2(p[1], p[0]));

                            boolean added = false;
                            for (int pIdx = 0; pIdx < nPeaks-1; pIdx++) {

                                if (wrap_PI(a-ap[pIdx])>0 && wrap_PI(a-ap[pIdx+1])<=0 && idxMapLocal[rIdx][x+d*y]==0) {

                                    idxMapLocal[rIdx][x+d*y] = nPeaks+1+pIdx;
                                    regSzsLocal[rIdx][nPeaks+1+pIdx]++;
                                    kernelLocal[rIdx][x+d*y] = -1;
                                    sumOFF++;
                                    added = true;
                                    //break;
                                }

                            }
                            if (!added) {
                                idxMapLocal[rIdx][x+d*y] = 2*nPeaks;
                                regSzsLocal[rIdx][2*nPeaks]++;
                                kernelLocal[rIdx][x+d*y] = -1;
                                sumOFF++;
                            }

                        }

                    }
                    else if ( d2 <= Rinner*Rinner ) {
                        idxMapLocal[rIdx][x+d*y] = 0;
                        regSzsLocal[rIdx][0]++;
                        kernelLocal[rIdx][x+d*y] = +1;
                        sumON++;
                    }
                    else {
                        idxMapLocal[rIdx][x+d*y] = -1; // nothing
                        kernelLocal[rIdx][x+d*y] = 0;
                    }
                }
            }

            System.out.println("\nrot. "+rIdx);
            for (int i = 0; i<2*nPeaks+1; i++) {
                    System.out.print("["+i+" -> "+regSzsLocal[rIdx][i]+"],");
            }
            System.out.println();

            // if it happens that it's too small region in one of the rotations, don't add the whole rotation
            for (int k = 0; k < 2*nPeaks+1; k++) {

                if (regSzsLocal[rIdx][k] <= 3) { // allow second sum to be 0 elements (overlaping with first)  && k!=2 && k!=4

                    System.out.println("some reg was small "+k+"th sum   : "+regSzsLocal[rIdx][k]+" elements");
                    return false;
                }

            }

            // normalize kernel for this rotation
            for (int k = 0; k<d*d; k++) {
                if (kernelLocal[rIdx][k]>0) {
                    kernelLocal[rIdx][k] /= sumON;
                }
                else if (kernelLocal[rIdx][k]<0) {
                    kernelLocal[rIdx][k] /= sumOFF;
                }
            }

        }

        index_map[0] = idxMapLocal;
        region_size[0] = regSzsLocal;
        kernels_rot[0] = kernelLocal;

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

//    public ImageStack plot()
//    {
//
//        ImageStack is = new ImageStack(d, d);
//
//        ImageProcessor ip = new FloatProcessor(d, d, regionIdxMap.get(kernelIdx*nRot));
//        is.addSlice("kernel,"+names.get(kernelIdx), ip);
//
//        return is;
//
//    }


    public ImageStack plotTemplates()
    {

        ImageStack is = new ImageStack(d, d);

        for (int l = 0; l<angles.size(); l++) {

            ImageProcessor ip = new FloatProcessor(d, d, regionIdxMap.get(l*nRot));
            is.addSlice("template,"+names.get(l), ip);

        }

        return is;

    }

    public ImageStack plotTemplatesAll()
    {

        ImageStack is = new ImageStack(d, d);

        for (int l = 0; l<regionIdxMap.size(); l++) {

            ImageProcessor ip = new FloatProcessor(d, d, regionIdxMap.get(l));
            is.addSlice("template,"+names.get(l/nRot), ip);

        }

        return is;

    }

    public ImageStack plotKernels()
    {

        ImageStack is = new ImageStack(d, d);

        for (int l = 0; l<angles.size(); l++) {

            ImageProcessor ip = new FloatProcessor(d, d, kernels.get(l*nRot));
            is.addSlice("kernel,"+names.get(l), ip);

        }

        return is;

    }

    public ImageStack plotKernel(int configurationIdx)
    {

        ImageStack is = new ImageStack(d, d);


        ImageProcessor ip = new FloatProcessor(d, d, kernels.get(configurationIdx));
        is.addSlice("kernel,"+names.get(configurationIdx/nRot), ip);

//        for (int l = 0; l<angles.size(); l++) {
//        }

        return is;

    }

    public ImageProcessor fit(FloatProcessor inip){

        int W = inip.getWidth();
        int H = inip.getHeight();

        ImageProcessor ipIdxs = new FloatProcessor(inip.getWidth(), inip.getHeight());
        ImageProcessor ipScos = inip.duplicate();
        convolver.convolveFloat(ipScos, kernels.get(0), d, d);//new FloatProcessor(inip.getWidth(), inip.getHeight());

       // new ImagePlus("after.convolving", ipScos).show();

        // find the best fit configuration/rotation for this feature
        // number of convolutions
        for (int convIdx = 1; convIdx<kernels.size(); convIdx++) {

//            System.out.println("convolving "+convIdx+" ...");

            ImageProcessor ipconv = inip.duplicate();
            convolver.convolveFloat(ipconv, kernels.get(convIdx), d, d);
            // check if there was max
            for (int l = 0; l<H*W; l++) {
                if(ipconv.getf(l)>ipScos.getf(l)) {
                    ipScos.setf(l, ipconv.getf(l));
                    ipIdxs.setf(l, convIdx);
                }
            }
        }

        return  ipIdxs;

    }

}