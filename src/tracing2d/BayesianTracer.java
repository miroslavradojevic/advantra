package tracing2d;

/**
 * Created by miroslav on 16-12-14.
 */
public class BayesianTracer {

    public static float[]       radiuses;
    public static int[]         nsamp2;
    public static float[][]     tt; // templates
    public static float[]       tta;

    public static float cross_sampling = .5f;

    public static void init(int rmin, int rmax) {

        for (int i = rmin; i <= rmax; i++) {

            radiuses[i] = i;
            nsamp2[i] = Math.round((3 * radiuses[i]) / cross_sampling);

            tt[i] = new float[2*nsamp2[i]+1];
            tta[i] = 0;

            for (int j = 0; j < tt[i].length; j++) {
                tt[i][j] = ;
            }

        }



    }



}
