/*
 TODO
 - input nwords comma separated, take all the values and calculate centroids for each nwords without recomputing the features
 - k means in parallel, domestic production
 - export to the output directory with name defined by it's parameters, for example NMACHINE.D.nw.op.on_500_60_50_50
 - output directory has c.txt, legend.txt and feats.csv
 - save the tagmap for each mXY.tif
 - calculate histograms - those would be the final feature to submit to train with
 */
package featureextraction;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import java.awt.Point;
import java.awt.Rectangle;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import mpi.cbg.fly.Feature;
import mpi.cbg.fly.Filter;
import mpi.cbg.fly.FloatArray2D;
import mpi.cbg.fly.FloatArray2DSIFT;
import mpi.cbg.fly.FloatArray2DScaleOctave;
import mpi.cbg.fly.ImageArrayConverter;
import weka.clusterers.SimpleKMeans;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.EuclideanDistance;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.ArffSaver;

/**
 *
 * @author Gadea
 */
public class FE implements PlugIn {

    private ImagePlus inimg;
    public byte[] inimgarray;

    private String inputDirectory;
    private String outputDirectory;

    private int D; //dimension of the square (size of ann edge)
    private double percentOver; // overlap percentage to make a grid
    private int[] nwords; // numberS of words for Bag of Words algorithm (K-means + weka)
    private double percentPOS; //overlap percentage to comparer the patches to get the positive patches
    private boolean save = false;

    //parameters to calculate SIFT features
    // steps
    private static int steps = 3;
    // initial sigma
    private static float initial_sigma = 1.6f;
    // feature descriptor size
    private static int fdsize = 4;
    // feature descriptor orientation bins
    private static int fdbins = 8;
    // size restrictions for scale octaves, use octaves < max_size and > min_size only
    private static int min_size = 64;
    private static int max_size = 1024;
    private static boolean upscale = false;
    private static float scale = 1.0f;

    //------------------
    @Override
    public void run(String arg) {

        inputDirectory = Prefs.get("neuronmachine.inputDirectory", System.getProperty("user.home"));
        D = (int) Prefs.get("neuronmachine.D", 500);
        percentOver = (double) Prefs.get("neuronmachine.percentOver", 50);
        percentPOS = (double) Prefs.get("neuronmachine.percentPOS", 50);
        String words = Prefs.get("neuronmachine.words", "20,40,60,80");
//        save = (boolean) Prefs.get("neuronmachine.save", false);

        GenericDialog gdG = new GenericDialog("to_make_a_grid");
        gdG.addStringField("input_directory", inputDirectory, 80);
        gdG.addNumericField("square_size", D, 0);
        gdG.addNumericField("percentage of overlap", percentOver, 0);
        gdG.addNumericField("percentage to be Positive", percentPOS, 0);
        gdG.addStringField("number of words", "20", 0);
//        gdG.addCheckbox("save the patches", save);

        gdG.showDialog();
        if (gdG.wasCanceled()) {
            return;
        }
        inputDirectory = gdG.getNextString();
        Prefs.set("neuronmachine.inputDirectory", inputDirectory);
        D = (int) gdG.getNextNumber();
        Prefs.set("neuronmachine.D", D);
        percentOver = (double) gdG.getNextNumber();
        Prefs.set("neuronmachine.percentOver", percentOver);
        percentPOS = (double) gdG.getNextNumber();
        Prefs.set("neuronmachine.percentPOS", percentPOS);
        words = gdG.getNextString();
        Prefs.set("neuronmachine.nwords", words);
//        save = (boolean) gdG.getNextBoolean();
//        Prefs.set("neuronmachine.save", save);
        String[] nw = words.split(",");
        nwords = new int[nw.length];
        for (int i = 0; i < nw.length; i++) {
            nwords[i] = Integer.parseInt(nw[i]);
        }

        if (inputDirectory.isEmpty()) {
            return;
        }

        outputDirectory = inputDirectory + File.separator + "SIFT.D.og.op.nw." + D + "_" + percentOver + "_" + percentPOS + "_";

        int overlapPixels = (int) percentPOS * D * D / 100; //minimum number of pixels are neccessary to consider a patch like positive (neuron)

        File dir = new File(inputDirectory);
        File[] imgFiles = listFilesEndingWith(dir, ".tif");
        File[] logFiles = listFilesEndingWith(dir, ".log");

        ArrayList<Vector<Feature>> siftFeatures = new ArrayList<Vector<Feature>>();//it stores all sift-features
        ArrayList<Integer> numSiftFeat = new ArrayList<Integer>();//it stores the number of sift descriptors for patch
        ArrayList<int[]> locationsxy = new ArrayList<int[]>();
        ArrayList<String> mosaicFileName = new ArrayList<String>();
        ArrayList<Integer> classPatches = new ArrayList<Integer>();//0 if it's no neuron and 1 it's neuron

        int total_desriptors = 0;
        for (int i = 0; i < imgFiles.length; i++) {
            long startTime = System.currentTimeMillis();
            //read the annot file and extract the neuron of gold-standard
            ArrayList<String> type_annot = new ArrayList<String>();  // mosaic name
            ArrayList<Point> points_annot = new ArrayList<Point>(); //it stores the upper-left corner to generate a rectangle
            try {
                FileReader fr;
                BufferedReader br;
                String line;
                fr = new FileReader(logFiles[i].getPath());
                br = new BufferedReader(fr);

                br.readLine(); //to read the comment
                while ((line = br.readLine()) != null) {
                    String[] aux = line.split("\t", 5);//[0] TYPE, [1] x, [2] y, [3] D, [4] D
                    type_annot.add(aux[0]);
                    points_annot.add(new Point(Integer.parseInt(aux[1]), Integer.parseInt(aux[2])));
                }
            } catch (IOException ex) {
                System.out.println("ERROR: " + ex.getMessage());
            }
            //neurons[] stores the points which represent the neurons of gold-standard
            Rectangle[] neurons = new Rectangle[points_annot.size()];

            for (int j = 0; j < points_annot.size(); j++) {
                String actual_class = type_annot.get(j);
                if (actual_class.contains("NEU")) {
                    neurons[j] = new Rectangle((int) points_annot.get(j).getX(), (int) points_annot.get(j).getY(), D, D);
                }
            }

            //for each mosaic
            File imgF = imgFiles[i];
            inimg = new ImagePlus(imgF.getAbsolutePath());
            inimg.show();
            int W = inimg.getWidth();
            int H = inimg.getHeight();

            inimgarray = (byte[]) inimg.getProcessor().getPixels();

            //get a tag map (binary map) where the white pixels depict the upper-left corner of the patches which shared at least the overlap percentage of their area with a neuron marked by the expert.
            ImagePlus impMap = getTagMap(neurons, inimg, overlapPixels);
//            impMap.show();

            ImageProcessor ipMap = impMap.getProcessor();
            byte[] pixelsMap = (byte[]) ipMap.getPixels();

            String path = inputDirectory + "\\tagMap\\";
            File fpath = new File(path);
            if (!fpath.exists()) {
                fpath.mkdirs();
            }
            IJ.save(impMap, path + impMap.getShortTitle() + ".tiff");

            Overlay ov = new Overlay();
            FileSaver fs;

            int margin = D / 2; //the margin

            double auxPercentage = 100 - percentOver;
            double step = auxPercentage * D / 100;

            int count = 0;
            int auxnn = 0;
            Rectangle rec = new Rectangle();

            // -------------------- seq  --------------------
            ArrayList<int[]> locXY = new ArrayList<int[]>();
//            ArrayList<Vector<Feature>> vv=new ArrayList<Vector<Feature>>();
//            int[] auxXY = new int[2];
//            int avgmillsec=0;
//            int avgmillseccounter=0;
//            byte[] patchyeah = new byte[D*D];
//            int countrubbish = 0;
//            long t11=System.currentTimeMillis();

            for (int x = margin; x < W - margin - D; x += step) {
                for (int y = margin; y < H - margin - D; y += step) {
//                    countrubbish++;
//                    if (countrubbish>5) return;
//                    auxXY[0] = x;
//                    auxXY[1] = y;
                    locXY.add(new int[]{x, y});
//                    long t1=System.currentTimeMillis();
//                    Vector<Feature> feat = ParallelSift.calculateSIFT(inimg, x, y, D);
//                    IJ.log(x+", "+y+" ");
//                    Vector<Feature> feat = ParallelSift.calculateSIFT(inimgarray, W, x, y, D, patchyeah); 
//                    vv.add(feat);
//                    long t2= System.currentTimeMillis();
//                    avgmillseccounter++;
//                    avgmillsec+=t2-t1;
//                    IJ.log("took "+ ((t2-t1)/1000f) +" s.");
                }
            }
//            long t22=System.currentTimeMillis();
//            IJ.log("sequential took "+(t22-t11)/1000f);    
//            IJ.log("avgtime was "+(avgmillsec/avgmillseccounter)/1000f);

            int CPU_NR = Runtime.getRuntime().availableProcessors();

            IJ.log("calculating SIFT features... \n" + imgFiles[i]);
            long t1 = System.currentTimeMillis();
            ParallelSift.load(inimg, locXY, pixelsMap, D);
            ParallelSift pSjobs[] = new ParallelSift[CPU_NR];

            for (int iJob = 0; iJob < pSjobs.length; iJob++) {
                pSjobs[iJob] = new ParallelSift(iJob * locXY.size() / CPU_NR, (iJob + 1) * locXY.size() / CPU_NR);
//                IJ.log("first argument "+iJob*locXY.size()/CPU_NR+" second argument "+(iJob+1)*locXY.size()/CPU_NR+" size locXY "+locXY.size());
                pSjobs[iJob].start();
            }
            for (int iJob = 0; iJob < pSjobs.length; iJob++) {
                try {
                    pSjobs[iJob].join();
                } catch (InterruptedException ex) {
                    ex.printStackTrace();
                }
            }

            long t2 = System.currentTimeMillis();
            IJ.log("" + (t2 - t1) / 1000f + " s.");

            // concatenate the data per image
            for (int j = 0; j < ParallelSift.siftFeatures.size(); j++) {
                siftFeatures.add((Vector<Feature>) ParallelSift.siftFeatures.get(j).clone());
                total_desriptors += ParallelSift.siftFeatures.get(j).size();
                numSiftFeat.add(ParallelSift.siftFeatures.get(j).size());
                locationsxy.add(ParallelSift.locXY.get(j).clone());
                mosaicFileName.add(imgFiles[i].getName());
                classPatches.add(ParallelSift.classPatches.get(j));
            }

        }

        IJ.log(IJ.d2s(total_desriptors / 1000f, 1) + "k desriptors to cluster");
        IJ.log("Kmeans centroids (weka)...");
        long startTime = System.currentTimeMillis();
        Instances ins_SIFT = getInstancesClustering(siftFeatures);
//        IJ.log("done");
        for (int i = 0; i < nwords.length; i++) {
            //all sift features have been calculated for all patches
            //obtain the histogram per each patch using Bag of Words algorithm
            //to do this, k-means algorithm on weka is used
            //k-means
            String auxOutput = outputDirectory;
            auxOutput += nwords[i];
            IJ.log("clustering... " + nwords[i]);
            float[][] centroids = wekaSKM(ins_SIFT, nwords[i]); // Ncentroids x 128
            exportCentroids(centroids, auxOutput, "c");

            //get the histogram per each patch,each patch is expressed with historgram wrt number of words
            IJ.log("getting the histogram per each patch");

            int[][] feat = getHistograms(centroids, siftFeatures); // Npatches x nwords
            IJ.log("done.");

            // export freatures (here those are histograms per centroids for each patch)
            exportFeaturesFull(feat, numSiftFeat, locationsxy, mosaicFileName, classPatches, auxOutput, "feat");

            // export legend

        }

        long endTime = System.currentTimeMillis();
        IJ.log("step machine learning took" + " " + ((endTime - startTime) / 60000f) + " min.");
    }

    protected void createTxtFile() {
        String path = "";
        File TXTfile = new File(path);

    }

    //it isn't used here. create a xlsFile with the results
//    protected void createXlsFile(ArrayList<String> namePatches, ArrayList<Integer[]> histoPatches, ArrayList<Integer> classPatches, int nwords) throws FileNotFoundException, IOException {
//        String path = inputDirectory + "\\sift_histograms.xls";
//        File XLSfile = new File(path);
//        //if there is a file with the same name, it's deleted
//        if (XLSfile.exists()) {
//            XLSfile.delete();
//        }
//        XLSfile.createNewFile();
//        //to create a book in an excel-file
//        Workbook book = new HSSFWorkbook();
//        //to start the connection with the xls-file
//        FileOutputStream file = new FileOutputStream(XLSfile);
//        try {
//            ArrayList<String> listHeader = new ArrayList<String>();
//            for (int i = 0; i < nwords; i++) {
//                listHeader.add("histo[" + i + "]");
//            }
//            listHeader.add("class");
//            listHeader.add("namePatch");
//            //write the data
//            Sheet sheet = book.createSheet("Sheet" + 1);
//            int n = 0;
//
//            Row row = sheet.createRow(n);
//            for (int j = 0; j < listHeader.size(); j++) {
//                Cell cell = row.createCell(j);
//                cell.setCellValue(listHeader.get(j));
//            }
//            for (int i = 0; i < namePatches.size(); i++) {
//                n++;
//                row = sheet.createRow(n);
//                for (int j = 0; j < nwords; j++) {
//                    Cell cell = row.createCell(j);
//                    cell.setCellValue(histoPatches.get(i)[j]);
//                }
//                Cell cell = row.createCell(nwords);
//                cell.setCellValue(classPatches.get(i));
//                cell = row.createCell(nwords + 1);
//                cell.setCellValue(namePatches.get(i));
//            }
//
//        } finally {
//            //to write the new data in the book
//            book.write(file);
//            //to close the connection with the file
//            file.close();
//        }
//    }

    //get the histograms per each patch from the classification of the clusters
    private int[][] getHistograms(float[][] centroids, ArrayList<Vector<Feature>> sift_feat) {
        // centroids: nwords x 128
        int[][] feat = new int[sift_feat.size()][centroids.length];

        // todo: this should be parallel
        for (int i = 0; i < sift_feat.size(); i++) { // through patches
            for (int j = 0; j < sift_feat.get(i).size(); j++) { // through keypoint locs of the patch

                float min_dist = Float.POSITIVE_INFINITY;

                int min_idx = -1;

                // find the closest centroid index
                for (int k = 0; k <centroids.length; k++) {

                    float dist = 0; // euclidean distance used (squared eucl.)
                    for (int l = 0; l < centroids[k][l]; l++) {
                        dist += Math.pow(centroids[k][l]-sift_feat.get(i).get(j).descriptor[l], 2);
                    }
                    if (dist<min_dist) {
                        min_dist = dist;
                        min_idx = k;
                    }

                }

                feat[i][min_idx]++;

            }
        }

        return feat;
    }

    //cluster the data with k-means and weka
    protected float[][] wekaSKM(Instances ins_data, int nwords) {
        long startTime = System.currentTimeMillis();
//        ArrayList<Integer> listclusters = new ArrayList<Integer>();
        //the instance has 128 attributes
        //this filter removes the attribute class
        SimpleKMeans skm = new SimpleKMeans();
        skm.setSeed(42);//this seed is by default

        try {
            skm.setNumClusters(nwords);
            skm.buildClusterer(ins_data);
        } catch (Exception ex) {
            Logger.getLogger(FE.class.getName()).log(Level.SEVERE, null, ex);
        }

        Instances insCentroids = skm.getClusterCentroids();
        //save the centroids in an arff-file. To test later.
//            extractArffClustering(insCentroids, "centroids"); 

        //float[][] format
        float[][] cent = new float[insCentroids.size()][128];
        for (int i = 0; i < insCentroids.size(); i++) {
            for (int j = 0; j < insCentroids.get(i).numValues(); j++) {
                cent[i][j] = (float) insCentroids.get(i).value(j);
            }
        }

//        try {
//            skm.setNumClusters(nwords);
//            skm.setSeed(42);//this seed is by default
//            skm.buildClusterer(ins_data);
//
//            //save a model-file. The file is opened with Weka 3.7 (interface). To test later.
////            weka.core.SerializationHelper.write(inputDirectory + "\\_model_SKM" + ".model", skm);
//
//            //save the centroids in an arff-file. To test later.
//            Instances insCentroids = skm.getClusterCentroids();
//            extractArffClustering(insCentroids, "centroids");
//            
//            //float[][] format
//            float[][] cent= new float[insCentroids.size()][128];
//            for (int i = 0; i < insCentroids.size(); i++) {
//                for (int j = 0; j < insCentroids.get(i).numValues(); j++) {
//                    cent[i][j]=(float)insCentroids.get(i).value(j);
//                }
//            }
//            
//            exportCentroids(cent,inputDirectory,"c");
        //.csv
//            extractCentroidsToCSV(insCentroids, "centroids");
//            listclusters = calculateDistanceToClosestCentroid(insCentroids, ins_data, skm);
        //using the clusters assigned by default 
//            for (int i = 0; i < ins_data.size(); i++) {
//                int cluster = skm.clusterInstance(ins_data.get(i));
//                listclusters.add(cluster);
//            }
//        } catch (Exception ex) {
//            Logger.getLogger(FE.class.getName()).log(Level.SEVERE, null, ex);
//        }
//        IJ.log("the model, calculates with SKM(for BoW), has been saved in " + inputDirectory + "\\_model_SKM" + ".model");
        long endTime = System.currentTimeMillis();
        long totalTime = endTime - startTime;
        double auxTime = totalTime / 60000;
        IJ.log("clustering took " + auxTime + " min.");
        return cent;
    }

    //obtain the instances to work with weka
    protected Instances getInstancesClustering(ArrayList<Vector<Feature>> siftFeatures) {
//        long startTime = System.currentTimeMillis();
        int ncolumns = siftFeatures.get(0).get(0).descriptor.length;//number of features

        //create a feature vector and the attributes to the feature vector
        FastVector fvWekaAttributes = new FastVector(ncolumns);
        for (int i = 0; i < ncolumns; i++) {
            fvWekaAttributes.addElement(new Attribute(IJ.d2s(i, 0)));
        }
        Instances ins_data = new Instances("Rel", fvWekaAttributes, ncolumns);
        // We create the instances and add them to the set of instances
        for (int i = 0; i < siftFeatures.size(); i++) {
            for (int j = 0; j < siftFeatures.get(i).size(); j++) {
                Instance ins = new DenseInstance(ncolumns);
                for (int k = 0; k < siftFeatures.get(i).get(j).descriptor.length; k++) {
                    ins.setValue((Attribute) fvWekaAttributes.elementAt(k), siftFeatures.get(i).get(j).descriptor[k]);
                }
                ins_data.add(ins);
            }
        }
        return ins_data;
    }

    //to obtain which is the closest cluster per each instance of data. 
    protected ArrayList<Integer> calculateDistanceToClosestCentroid(Instances clusterCentroid, Instances data, SimpleKMeans skm) {
//        SimpleKMeans skm = new SimpleKMeans();
        EuclideanDistance Dist = (EuclideanDistance) skm.getDistanceFunction();
        ArrayList<Integer> clustersInst = new ArrayList<Integer>();//store the number of cluster for each instance
        int numCluster = -1;

        for (int iI = 0; iI < data.numInstances(); iI++) {
            double minDistance = Dist.distance(clusterCentroid.get(0), data.instance(iI));
            numCluster = 0;
            for (int iC = 1; iC < clusterCentroid.size(); iC++) {
                Double dist = Dist.distance(clusterCentroid.get(iC), data.instance(iI));
                if (dist < minDistance) {
                    minDistance = dist;
                    numCluster = iC;
                }
            }
            clustersInst.add(numCluster);
        }
        return clustersInst; //return an array per each instance with the number of cluster which is the closest.
    }

    //get the sift features per each patch
    protected Vector<Feature> calculateSIFT(ImagePlus imp) {
        ImageProcessor ip1 = imp.getProcessor().convertToFloat();

        Vector< Feature> fs1;
        FloatArray2DSIFT sift = new FloatArray2DSIFT(fdsize, fdbins);
        FloatArray2D fa = ImageArrayConverter.ImageToFloatArray2D(ip1);
        Filter.enhance(fa, 1.0f);

        if (upscale) {
            FloatArray2D fat = new FloatArray2D(fa.width * 2 - 1, fa.height * 2 - 1);
            FloatArray2DScaleOctave.upsample(fa, fat);
            fa = fat;
            fa = Filter.computeGaussianFastMirror(fa, (float) Math.sqrt(initial_sigma * initial_sigma - 1.0));
        } else {
            fa = Filter.computeGaussianFastMirror(fa, (float) Math.sqrt(initial_sigma * initial_sigma - 0.25));
        }
        sift.init(fa, steps, initial_sigma, min_size, max_size);
        fs1 = sift.run(max_size);
        Collections.sort(fs1);
        return fs1;
    }

    public String getExtension(String filePath) {
        String extension = "";
        int dotIdx = filePath.lastIndexOf('.');
        if (dotIdx > 0) {
            extension = filePath.substring(dotIdx + 1);
        }
        return extension;
    }

    public String removeExtension(String filePath) {
        String extension = "";
        int dotIdx = filePath.lastIndexOf('.');
        if (dotIdx > 0) {
            extension = filePath.substring(0, dotIdx);
        }
        return extension;
    }

    public File[] listFilesEndingWith(File dir, String suffix) {
        final String sfx = suffix;
        File[] tif_train = dir.listFiles(
                new FilenameFilter() {
                    public boolean accept(File dir, String name) {
                        return name.toLowerCase().endsWith(sfx);
                    }
                }
        );
        return tif_train;
    }

    private double overlap(Rectangle rect1, Rectangle rect2) {

        double x11 = rect1.getX();//left -up
        double y11 = rect1.getY();//left -up

        double x12 = rect1.getX() + rect1.getWidth();//right -up

        double y12 = rect1.getY() + rect1.getHeight();//left-down

        double x21 = rect2.getX();//left -up
        double y21 = rect2.getY();//left -up

        double x22 = rect2.getX() + rect2.getWidth();//right -up

        double y22 = rect2.getY() + rect2.getHeight();//left-down

        double xDiff;
        if (x11 < x21) {
            xDiff = x12 - x21;
        } else {
            xDiff = x22 - x11;
        }

        double yDiff;
        if (y11 < y21) {
            yDiff = y12 - y21;
        } else {
            yDiff = y22 - y11;
        }
        xDiff = (xDiff < 0) ? 0 : xDiff;
        yDiff = (yDiff < 0) ? 0 : yDiff;

        return xDiff * yDiff;
    }

    public void extractCentroidsToCSV(Instances instances, String name) {
        File file = new File(inputDirectory + File.separator + name + ".csv");
        FileWriter fw = null;
        try {
            fw = new FileWriter(file.getAbsolutePath());
            String header = "d[0]";
            for (int i = 1; i < instances.get(0).numAttributes(); i++) {
                header += ",d[" + i + "]";
            }
            fw.append(header);
            fw.append("\n");
            String centroid;
            for (int i = 0; i < instances.size(); i++) {
                centroid = String.valueOf(instances.get(i).value(0));
                for (int j = 1; j < instances.get(i).numAttributes(); j++) {
                    centroid += "," + instances.get(i).value(j);
                }
                fw.append(centroid);
                fw.append("\n");
            }
        } catch (Exception e) {
            e.getMessage();
        } finally {
            try {
                fw.flush();
                fw.close();
            } catch (IOException ex) {
                Logger.getLogger(FE.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
    }

    //to extract the instances to arff-file.
    public void extractArffClustering(Instances instances, String name) {
        ArffSaver saver = new ArffSaver();
        saver.setInstances(instances);
        try {
            saver.setFile(new File(inputDirectory + File.separator + name + ".arff"));
            saver.writeBatch();
        } catch (Exception ex) {
            IJ.log(ex.getMessage());
        }
    }

    protected ImagePlus getTagMap(Rectangle[] neurons, ImagePlus inimg, double overlapPixels) {
        ImagePlus impMap = IJ.createImage(inimg.getShortTitle() + "_map", "8-bit black", inimg.getWidth(), inimg.getHeight(), 1);
        ImageProcessor ipMap = impMap.getProcessor();
        byte[] pixelsMap = (byte[]) ipMap.getPixels();

        Overlay ovMap2 = new Overlay();
        Rectangle aux = new Rectangle();
        for (Rectangle neuron : neurons) {
            Rectangle recN = neuron;
            Roi rec_roi = new Roi(recN);
            ovMap2.add(rec_roi);
            impMap.setOverlay(ovMap2);

            int x1 = (int) recN.getX() - D; //upper-left corner of the neuron
            int y1 = (int) recN.getY() - D;
            int x2 = (int) recN.getX() + D; //upper-left corner of the neuron
            int y2 = (int) recN.getY() + D;

            for (int i = x1; i <= x2; i++) {
                for (int j = y1; j <= y2; j++) {
                    if ((i >= 0) && (j >= 0) && (i < inimg.getWidth()) && (j < inimg.getHeight())) {
                        aux.setRect(i, j, D, D);
                        double overAux = overlap(aux, neuron);
                        if (overAux >= overlapPixels) {
                            int pos = j * inimg.getWidth() + i;
                            pixelsMap[pos] = (byte) 255;
                        }
                    }
                }
            }

        }
        ipMap.setPixels(pixelsMap);
        impMap.setProcessor(ipMap);
        return impMap;
    }

    //Gadea's version
    protected ImagePlus getTagMap_1(Rectangle[] neurons, ImagePlus inimg) {
        ImagePlus impMap = IJ.createImage(inimg.getShortTitle() + "_map", "8-bit black", inimg.getWidth(), inimg.getHeight(), 1);
        ImageProcessor ipMap = impMap.getProcessor();
        byte[] pixelsMap = (byte[]) ipMap.getPixels();
        Overlay ovMap2 = new Overlay();

        for (Rectangle neuron : neurons) {
            Rectangle recN = neuron;
            Roi rec_roi = new Roi(recN);
            ovMap2.add(rec_roi);
            impMap.setOverlay(ovMap2);

            double auxPercentPOS = percentPOS * D / 100;
            int d = (int) (D - auxPercentPOS);

            int r2 = (D / 2 + d) * (D / 2 + d) + D / 2 * D / 2; //squared ratio of circunference

            int x = (int) recN.getX(); //upper-left corner of the neuron
            int y = (int) recN.getY();

            int centerX = (int) (x + D / 2); //centroid of the neuron
            int centerY = (int) (y + D / 2);

            //upper-left of circular sector
            int cX = centerX;
            int cY = centerY;
            for (int i = x - d; i < x; i++) {
                for (int j = y - d; j < y; j++) {
                    if (((i - cX) * (i - cX) + (j - cY) * (j - cY)) <= r2) {
                        if ((i >= 0) && (i < impMap.getWidth()) && (j >= 0) && (j < impMap.getHeight())) {
                            int offset = j * impMap.getWidth();
                            int pos = offset + i;
                            pixelsMap[pos] = (byte) 255;
                            impMap.updateAndDraw();
                        }
                    }
                }
            }
            //upper-rigth of circular sector
            cX = centerX - 2 * d;
            cY = centerY;
            for (int i = x; i < x + d; i++) {
                for (int j = y - d; j < y; j++) {
                    if (((i - cX) * (i - cX) + (j - cY) * (j - cY)) <= r2) {
                        if ((i >= 0) && (i < impMap.getWidth()) && (j >= 0) && (j < impMap.getHeight())) {
                            int offset = j * impMap.getWidth();
                            int pos = offset + i;
                            pixelsMap[pos] = (byte) 255;
                            impMap.updateAndDraw();
                        }
                    }
                }
            }
            //bottom-right of circular sector
            cX = centerX - 2 * d;
            cY = centerY - 2 * d;
            for (int i = x; i < x + d; i++) {
                for (int j = y; j < y + d; j++) {
                    if (((i - cX) * (i - cX) + (j - cY) * (j - cY)) <= r2) {
                        if ((i >= 0) && (i < impMap.getWidth()) && (j >= 0) && (j < impMap.getHeight())) {
                            int offset = j * impMap.getWidth();
                            int pos = offset + i;
                            pixelsMap[pos] = (byte) 255;
                            impMap.updateAndDraw();
                        }
                    }
                }
            }
            //bottom-left of circular sector
            cX = centerX;
            cY = centerY - 2 * d;
            for (int i = x - d; i < x; i++) {
                for (int j = y; j < y + d; j++) {
                    if (((i - cX) * (i - cX) + (j - cY) * (j - cY)) <= r2) {
                        if ((i >= 0) && (i < impMap.getWidth()) && (j >= 0) && (j < impMap.getHeight())) {
                            int offset = j * impMap.getWidth();
                            int pos = offset + i;
                            pixelsMap[pos] = (byte) 255;
                            impMap.updateAndDraw();
                        }
                    }
                }
            }
            impMap.updateAndDraw();
        }
        return impMap;
    }

    public static void createAndCleanDir(String dirpath) {
        File f1 = new File(dirpath);
        if (!f1.exists()) {
            f1.mkdirs();
        } else { // clean it, reset with each exec
            File[] files = f1.listFiles();
            for (File file : files) {
                if (!file.delete()) {
                    System.out.println("Failed to delete " + file);
                }
            }
        }
//        IJ.log(dirpath);
    }

    private void exportCentroids(float[][] centroids, String outdir, String filename) {

        createDir(outdir);
        String outfile = outdir + File.separator + filename + ".txt";
        cleanfile(outfile);

        try {

            PrintWriter logWriter1 = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));

            int t1 = 0;

            for (int i = 0; i < centroids.length; i++) {
                if (centroids[i] != null) {
                    for (int j = 0; j < centroids[i].length; j++) {
                        logWriter1.print(centroids[i][j] + " ");
                    }
                    logWriter1.println("");
                }
            }

            logWriter1.close();

        } catch (IOException e) {
        }

    }

//    private void exportFeatures(int[][] feat, String outdir, String filename) {
//
//        createDir(outdir);
//        String outfile = outdir + File.separator + filename + ".csv";
//        cleanfile(outfile);
//
//        try {
//
//            PrintWriter logWriter1 = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));
//
//            int t1 = 0;
//
//            for (int i = 0; i < feat.length; i++) {
//                    for (int j = 0; j < feat[i].length; j++) {
//                        logWriter1.print(feat[i][j] + ",");
//                    }
//                    logWriter1.println("");
//            }
//
//            logWriter1.close();
//
//        } catch (IOException e) {
//        }
//
//        // export the legend at the same time with features
//
//    }

    private void exportFeaturesFull(
            int[][] feat,
            ArrayList<Integer> nkeypointlocs,
            ArrayList<int[]> locsXY,
            ArrayList<String> mosaicname,
            ArrayList<Integer> classpatch,
            String outdir, String filename) {

        createDir(outdir);
        String outfile = outdir + File.separator + filename + ".csv";
        cleanfile(outfile);

        try {

            PrintWriter logWriter1 = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));

            int t1 = 0;

            for (int i = 0; i < feat.length; i++) {

                for (int j = 0; j < feat[i].length; j++) logWriter1.print(feat[i][j] + ",");

                logWriter1.print(locsXY.get(i)[0]+","+locsXY.get(i)[1]+",");

                logWriter1.print(mosaicname.get(i)+",");

                logWriter1.println(Integer.toString(classpatch.get(i)));
            }

            logWriter1.close();

        } catch (IOException e) {}

        outfile = outdir + File.separator + filename + ".csv" + ".legend";
        cleanfile(outfile);

        try {

            PrintWriter logWriter1 = new PrintWriter(new BufferedWriter(new FileWriter(outfile, true)));

            logWriter1.println("# attributes, description of the dataset, each row has explanation for one of the features");

            int countfeat = 0;
            for (int j = 0; j < feat[j].length; j++) logWriter1.println(IJ.d2s(++countfeat,0)+ " " +"h"+IJ.d2s(countfeat,0));
            logWriter1.println(IJ.d2s(++countfeat, 0) + " " + "x");
            logWriter1.println(IJ.d2s(++countfeat, 0) + " " + "y");
            logWriter1.print(IJ.d2s(++countfeat, 0) + " " + "mosaic");
            logWriter1.println(IJ.d2s(++countfeat, 0) + " " + "classpatch");

            logWriter1.close();

        } catch (IOException e) {}

    }

    public static void createDir(String dirpath) {
        // create directory without cleaning it up
        File f1 = new File(dirpath);
        if (!f1.exists()) {
            f1.mkdirs();
        }
//        IJ.log("createDir " + dirpath);
    }

    public static void cleanfile(String filepath) {

        try {
            PrintWriter logWriter = new PrintWriter(filepath);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {
        }

    }

    public float[][] readCentroids(String pathToCentroidsCsv) {
        // 
        String csvFilePath = new File(pathToCentroidsCsv).getAbsolutePath();// path to csv file with centroids

        if (!(new File(csvFilePath).exists())) {
            IJ.log(csvFilePath + " does not exist!");
            return null;
        }

        IJ.log("reading...\t" + csvFilePath);

        ArrayList<float[]> centroids = new ArrayList<float[]>(); // 1x128 in each row (sift format)
        try { // scan the file
            FileInputStream fstream = new FileInputStream(csvFilePath);
            BufferedReader br = new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
            String read_line;
            while ((read_line = br.readLine()) != null) {
                if (read_line.isEmpty()) {
                    continue;
                }
                if (!read_line.trim().startsWith("#")) { // # are comments
                    String[] readLn = read_line.trim().replaceAll(",", ".").split("\\s+");
                    if (readLn.length != 128) {
                        continue; // skip the line that did not have enough values
                    }
                    float[] valsLn = new float[128]; // sift desriptor
                    for (int i = 0; i < 128; i++) {
                        valsLn[i] = Float.valueOf(readLn[i].trim()).floatValue();
                    }
                    centroids.add(valsLn);
                }

            }
            br.close();
            fstream.close();
        } catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
        }

        IJ.log(centroids.size() + " centroids read");

        // prepare ouput
        float[][] out = new float[centroids.size()][128];
        for (int i = 0; i < out.length; i++) {
            for (int j = 0; j < out[i].length; j++) {
                out[i][j] = centroids.get(i)[j];
            }
        }

        centroids.clear();
        return out;

    }

}