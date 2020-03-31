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
public class KmeansDemo implements PlugIn {

    public void run(String s) {
        IJ.log("*** KMeans demo ***");

        IJ.log("2d example");
        Random rand = new Random();
        int N = 50000;    // number of data vectors
        int M = 2;      // data vector dimensionality
        int K = 20;    // nr. clusters
        int MAX_ITER = 100;
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

        float[][] c = new float[K][data[0].length];

        IJ.log(K+"-means clustering...");
        long t1 = System.currentTimeMillis();
        kmeans(data, K, MAX_ITER, c);
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


//        if(true)return;
        //**************************************************************************************
        String feat_csv_path = Prefs.get("neuronmachine.feat_csv_path", "");
        GenericDialog gdG = new GenericDialog("Train");
        gdG.addStringField("input_directory", feat_csv_path, 80);
        gdG.showDialog();
        if (gdG.wasCanceled()) return;
        feat_csv_path = gdG.getNextString();
        Prefs.set("neuronmachine.feat_csv_path", feat_csv_path);

        File fl  = new File(feat_csv_path);
        if (!fl.exists()) return;

        ArrayList<double[]> feat = new ArrayList<double[]>();

        BufferedReader br = null;
        String line = "";
        String cvsSplitBy = ",";

        try {

            br = new BufferedReader(new FileReader(feat_csv_path));
            while ((line = br.readLine()) != null) {

                String[] datarow_string = line.split(cvsSplitBy);
                if (datarow_string.length<128) continue;

                double[] datarow = new double[128];
                for (int i = 0; i < datarow.length; i++)
                    datarow[i] = Double.valueOf(datarow_string[i].trim());

                feat.add(datarow);

            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) {
                try {
                    br.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        IJ.log("Done reading : " + feat.size() + " x " + feat.get(0).length);

        // convert to float[][] required by kmeans()
        K = 30;
        data = new float[feat.size()][feat.get(0).length];
        for (int i = 0; i < feat.size(); i++) {
            for (int j = 0; j < feat.get(i).length; j++) {
                data[i][j] = (float) feat.get(i)[j];
            }
        }
        c = new float[K][feat.get(0).length];

        IJ.log(K+"-means clustering...");
        t1 = System.currentTimeMillis();
        kmeans(data, K, MAX_ITER, c);
        t2 = System.currentTimeMillis();
        IJ.log("done. " + (t2 - t1) / 1000f + " s.");

    }

    private void kmeans(float[][] _data, int _K, int MAX_ITER, float[][] _c){ // _t are centroid assignments, _c are centroids

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
                try {jobs[iJ].join();} catch (InterruptedException ex) {ex.printStackTrace();}
            }

        }
        while (ThreadedKMeans.assignments_changed() && iter<MAX_ITER);

        if (iter==MAX_ITER) IJ.log("reached MAX_ITER="+MAX_ITER);

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