/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package kmeans;

import ij.IJ;
import ij.plugin.PlugIn;
import java.util.Arrays;
import java.util.Random;

/**
 *
 * @author Gadea
 */
public class Kmeans {

    private int N;    // number of data vectors
    private int M;      // data vector dimensionality
    private int K;    // nr. clusters
    private float[][] data;

    public Kmeans(int N, int M, int K) {
        this.N = N;
        this.M = M;
        this.K = K;
        this.data = new float[N][M];// generate data vectors
    }

    public float[][] run(float[][] data) {

        this.data = data.clone();// generate data vectors

        IJ.log("" + N + " " + M + "-dimensional vectors ");

        // clustering outputs:  - assignments int[] t, where t[i]\in[0,K)
        //                      - centroids float[][] c, where c[i][...] represents one centroid
        int[] t = new int[data.length];
        float[][] c = new float[K][data[0].length];

        IJ.log(K + "-means clustering...");
        long t1 = System.currentTimeMillis();
        kmeans(data, K, t, c);
        long t2 = System.currentTimeMillis();
        IJ.log("done. " + (t2 - t1) / 1000f + " s.");
        return c;
    }

    private void kmeans(float[][] _data, int _K, int[] _t, float[][] _c) { // _t are centroid assignments, _c are centroids

        // uses threaded kmenans clustering
        int CPU_NR = Runtime.getRuntime().availableProcessors();

        ThreadedKMeans.initialize(_data, _K);

        int iter = 0;
        do {

//            IJ.log("."); // "--- iter "+(++iter)
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

        } while (ThreadedKMeans.assignments_changed());

        // add the values to output
        for (int i = 0; i < ThreadedKMeans.t0.length; i++) {
            _t[i] = ThreadedKMeans.t0[i];
        }

        for (int i = 0; i < ThreadedKMeans.c0.length; i++) {
            for (int j = 0; j < ThreadedKMeans.c0[i].length; j++) {
                _c[i][j] = ThreadedKMeans.c0[i][j];
            }
        }

    }

}
