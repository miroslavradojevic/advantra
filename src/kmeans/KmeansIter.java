package kmeans;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 * Created by miroslav on 11-4-16.
 */
public class KmeansIter {

    private int N;    // number of data vectors
    private int M;      // data vector dimensionality
    private int K;    // nr. clusters
    private float[][] data;
    private int MAX_ITER = 100;

    public KmeansIter(int N, int M, int K, int MAX_ITER) {
        this.N = N;
        this.M = M;
        this.K = K;
        this.data = new float[N][M];// generate data vectors
        this.MAX_ITER = MAX_ITER;
    }

    public float[][] run(float[][] data) {
//        IJ.log("*** KMeans demo ***");
        this.data = data.clone();// generate data vectors
        IJ.log("" + N + " " + M + "-dimensional vectors ");

//        IJ.log("2d example");
        Random rand = new Random();
//        int N = 50000;    // number of data vectors
//        int M = 2;      // data vector dimensionality
//        int K = 20;    // nr. clusters
//        int MAX_ITER = 100;
//        float[][] data = new float[N][M];// generate data vectors
        float[] data_min = new float[M];
        Arrays.fill(data_min, Float.POSITIVE_INFINITY);
        float[] data_max = new float[M];
        Arrays.fill(data_max, Float.NEGATIVE_INFINITY);

        IJ.log("" + N + " " + M + "-dimensional vectors ");

        // clustering outputs:  - assignments int[] t, where t[i]\in[0,K)
        //                      - centroids float[][] c, where c[i][...] represents one centroid
        float[][] c = new float[K][data[0].length];

        IJ.log(K + "-means clustering...");
        long t1 = System.currentTimeMillis();
        kmeans(data, K, MAX_ITER, c);
        long t2 = System.currentTimeMillis();
        IJ.log("done. " + (t2 - t1) / 1000f + " s.");

        return c;

    }

    private void kmeans(float[][] _data, int _K, int MAX_ITER, float[][] _c) { // _t are centroid assignments, _c are centroids

        // uses ThreadedKMeans class for clustering
        int CPU_NR = Runtime.getRuntime().availableProcessors();

        ThreadedKMeans.initialize(_data, _K);

        int iter = 0;
        do {

//            IJ.log("iter="+iter);
            ++iter;

            ThreadedKMeans jobs[] = new ThreadedKMeans[CPU_NR];

            for (int iJ = 0; iJ < jobs.length; iJ++) {
                jobs[iJ] = new ThreadedKMeans(iJ * _data.length / CPU_NR, (iJ + 1) * _data.length / CPU_NR);
                jobs[iJ].start();
            }

            for (int iJ = 0; iJ < jobs.length; iJ++) {
                try {
                    jobs[iJ].join();
                } catch (InterruptedException ex) {
                    ex.printStackTrace();
                }
            }

        } while (ThreadedKMeans.assignments_changed() && iter < MAX_ITER);

        if (iter == MAX_ITER) {
            IJ.log("reached MAX_ITER=" + MAX_ITER);
        }

//        // add the values to output
//        for (int i = 0; i < ThreadedKMeans.t0.length; i++) {
//            _t[i] = ThreadedKMeans.t0[i];
//        }
        // c0 is given as output, c0<-c1 is assigned before that in assignments_changed()
        for (int i = 0; i < ThreadedKMeans.c0.length; i++) {
            for (int j = 0; j < ThreadedKMeans.c0[i].length; j++) {
                _c[i][j] = ThreadedKMeans.c0[i][j];
            }
        }

//        //*** DEBUG: export calculated centroids after number of iterations
//        String outfile = System.getProperty("user.home")+File.separator+"c.csv";
//        IJ.log(outfile);
//        String centroid_out = "";
//
//        try {
//            // output centroids
//            FileWriter writer = new FileWriter(outfile);
//            for (int i = 0; i < ThreadedKMeans.c0.length; i++) {
//                for (int j = 0; j < ThreadedKMeans.c0[i].length; j++) {
//                    centroid_out+=String.valueOf(ThreadedKMeans.c0[i][j])+((j<ThreadedKMeans.c0[i].length-1)?",":"\n");
//                }
//            }
//            writer.write(centroid_out);
//            writer.close();
//
//        }
//        catch(IOException e) {e.printStackTrace();}
//        //*** DEBUG
    }

}
