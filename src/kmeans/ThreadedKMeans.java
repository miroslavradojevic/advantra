package kmeans;

import ij.IJ;

import java.util.Arrays;
import java.util.Random;

/**
 * Created by miroslav on 11-4-16.
 */
public class ThreadedKMeans extends Thread {

    private int begN, endN;
    public static float[][] data; // data rows N * M
    public static float[][] data_range; // M * 2 (min--max)
    public static int K;

    public static float[][] c0; // centroids
    public static int[] t0; // clusters [0,K) values (previous)
    public static int[] t1; // clusters (current)

    // mutable objects shared by multiple threads, need to be modified using synchronized methods
    // http://www.programcreek.com/2014/02/how-to-make-a-method-thread-safe-in-java/
    public static float[][] c1;   // aux. centroids iterative avg.
    public static int[]     n1;  // aux. centroids count cluster size

    public ThreadedKMeans(int n0, int n1) {
        begN = n0;
        endN = n1;
    }

    public static void initialize(float[][] _data, int _K) {

//        iter = 0;

        data = _data;// try to keep the pointer only
        data_range = new float[_data[0].length][2];
        for (int i = 0; i < _data[0].length; i++) {
            data_range[i][0] = _data[0][i]; // min
            data_range[i][1] = _data[0][i]; // max
            for (int j = 1; j < _data.length; j++) {
                if (_data[j][i]<data_range[i][0]) data_range[i][0] = _data[j][i];
                if (_data[j][i]>data_range[i][1]) data_range[i][1] = _data[j][i];
            }
        }

        K = _K;
        c0 = new float[K][_data[0].length]; // centroids
        Random rnd = new Random();
        for (int i = 0; i < _data[0].length; i++) { // initial random centroids
            for (int j = 0; j < K; j++) {
                c0[j][i] = data_range[i][0] + rnd.nextFloat() * (data_range[i][1]-data_range[i][0]);
            }
        }

        t0 = null; // at the initialization so that first t0==t1 check is skipped
        t1 = new int[_data.length];

        // will need to be thread synchronized
        c1  = new float[K][_data[0].length]; // aux. centroid calculation storage
        n1 = new int[K]; // store number of elements in that cluster

    }

    public void run() {

        for (int lidx = begN; lidx < endN; lidx++) {

            float d2min = Float.MAX_VALUE;
            int   kmin = -1;

            for (int kidx = 0; kidx < K; kidx++) {

                float d2 = 0;
                for (int j = 0; j < data[0].length; j++) {
                    d2 += Math.pow(data[lidx][j]-c0[kidx][j],2);
                }

                if (d2<d2min) {
                    d2min = d2;
                    kmin = kidx;
                }

            }

            t1[lidx] = kmin;

            // use thread synchronization when calculating iterative mean in c1[kmin], n1[kmin]
            iterative_mean(kmin, data[lidx]); // update the centroids for the next iteration

        }

    }

    private static synchronized void iterative_mean(int _k, float[] _datarow){

        n1[_k]++;
        for (int i = 0; i < _datarow.length; i++) {
            c1[_k][i] = ((float)(n1[_k]-1)/n1[_k]) * c1[_k][i] + (1f/n1[_k]) * _datarow[i];
        }

    }

    public static boolean assignments_changed(){

        boolean tags_equal;

        if(t0==null) {
            t0 = new int[t1.length];
            tags_equal = false;
        }
        else{
            tags_equal = true;
            for (int i = 0; i < data.length; i++) {
                tags_equal = tags_equal && (t1[i]==t0[i]);
                if (!tags_equal) break;
            }
        }

        // t0 <- t1
        for (int i = 0; i < t0.length; i++) {
            t0[i] = t1[i];
        }

        // c0 <- c1, and reset n1, c1 to zero
        for (int i = 0; i < c0.length; i++) {
            n1[i] = 0;
            for (int j = 0; j < c0[i].length; j++) {
                c0[i][j] = c1[i][j];
                c1[i][j] = 0;
            }
//            IJ.log(Arrays.toString(c0[i]));
        }

        return !tags_equal;
    }

}
