package kmeans;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;

import java.awt.*;
import java.util.Arrays;
import java.util.Random;

/**
 * Created by miroslav on 11-4-16.
 */
public class KmeansDemo implements PlugIn {

    public void run(String s) {
        IJ.log("*** KMeans demo ***");

        IJ.log("2d example");
        Random rand = new Random();
        int N = 50000;    // number of data vectors
        int M = 2;      // data vector dimensionality
        int K = 20;    // nr. clusters
        float[][] data = new float[N][M];// generate data vectors
        float[] data_min = new float[M];
        Arrays.fill(data_min, Float.POSITIVE_INFINITY);
        float[] data_max = new float[M];
        Arrays.fill(data_max, Float.NEGATIVE_INFINITY);

        for (int k = 0; k < K; k++) {

            float[] mean = new float[M];
            float[] std  = new float[M];
            for (int j = 0; j < M; j++) {
                mean[j] = rand.nextFloat();
                std[j]  = 0.2f*rand.nextFloat();
            }

            for (int i = k * N / K; i < (k + 1) * N / K; i++) {
                for (int j = 0; j < M; j++) {
                    data[i][j] = (float) (mean[j] + std[j] * rand.nextGaussian());
                    if (data[i][j]>data_max[j]) data_max[j] = data[i][j];
                    if (data[i][j]<data_min[j]) data_min[j] = data[i][j];
                }
            }

        }

        IJ.log("" + N + " " + M + "-dimensional vectors ");

        // clustering outputs:  - assignments int[] t, where t[i]\in[0,K)
        //                      - centroids float[][] c, where c[i][...] represents one centroid

        int[] t = new int[data.length];
        float[][] c = new float[K][data[0].length];

        IJ.log(K+"-means clustering...");
        long t1 = System.currentTimeMillis();
        kmeans(data, K, t, c);
        long t2 = System.currentTimeMillis();
        IJ.log("done. " + (t2 - t1) / 1000f + " s.");

        // viz
        ImagePlus imp = new ImagePlus("", new ByteProcessor(512, 512));
        imp.show();
        Overlay ov = new Overlay();
        for (int i = 0; i < N; i++) {
            OvalRoi pr = new OvalRoi(
                    ((data[i][0]-data_min[0])/(data_max[0]-data_min[0]))*512f-.5f,
                    ((data[i][1]-data_min[1])/(data_max[1]-data_min[1]))*512f-.5f, 1f, 1f);
            pr.setFillColor(new Color(0,0,1f,.2f));
            ov.add(pr);
        }
        for (int i = 0; i < K; i++) {
            OvalRoi dummy = new OvalRoi(
                    ((c[i][0]-data_min[0])/(data_max[0]-data_min[0]))*512f-5f,
                    ((c[i][1]-data_min[1])/(data_max[1]-data_min[1]))*512f-5f, 10f, 10f);
//                    c[i][0]*512f-5f,
//                    c[i][1]*512f-5f, 10f, 10f);
            dummy.setFillColor(Color.RED);
            ov.add(dummy);
        }
        imp.setOverlay(ov);

    }

    private void kmeans(float[][] _data, int _K, int[] _t, float[][] _c){ // _t are centroid assignments, _c are centroids

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
                try {jobs[iJ].join();} catch (InterruptedException ex) {ex.printStackTrace();}
            }

        }
        while (ThreadedKMeans.assignments_changed());

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