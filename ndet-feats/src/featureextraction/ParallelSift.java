package featureextraction;

import ij.ImagePlus;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Vector;
import mpi.cbg.fly.Feature;
import mpi.cbg.fly.Filter;
import mpi.cbg.fly.FloatArray2D;
import mpi.cbg.fly.FloatArray2DSIFT;
import mpi.cbg.fly.FloatArray2DScaleOctave;
import mpi.cbg.fly.ImageArrayConverter;

/**
 *
 * @author Gadea
 */
public class ParallelSift extends Thread {

    private int begN, endN;

    public static ArrayList<Vector<Feature>> siftFeatures = new ArrayList<Vector<Feature>>();
    public static ImagePlus inimg;
    public static byte[] pixelsMap;
    public static byte[] inimgarray;
    public static ArrayList<int[]> locXY = new ArrayList<int[]>();
    public static int D;
    //    private static ArrayList<String> namePatches = new ArrayList<String>();
    public static ArrayList<Integer> classPatches = new ArrayList<Integer>();

    //parameters to calculate SIFT features
    // steps
    public static int steps;// = 3;
    // initial sigma
    public static float inSigma;// = 1.6f;
    // feature descriptor size
    public static int fdsize;// = 4;
    // feature descriptor orientation bins
    public static int fdbins;// = 8;
    // size restrictions for scale octaves, use octaves < max_size and > min_size only
    public static int min_size;// = 64;
    public static int max_size;// = 1024;
    public static boolean upscale;// = false;
    private static float scale = 1.0f;

    //------------------
    public ParallelSift(int n0, int n1) {
        begN = n0;
        endN = n1;
    }
    public static void load(ImagePlus inputImage, ArrayList<int[]> inputLocXY, byte[] inputMap, int inputD,int inpSteps,float inpSigma, int inpfdsize,int inpfdbins, int inpmin_size, int inpmax_size,boolean inpupscale) {
        siftFeatures.clear();
        locXY.clear();
        classPatches.clear();

        for (int i = 0; i < inputLocXY.size(); i++) {
            siftFeatures.add(new Vector<Feature>());
            locXY.add(inputLocXY.get(i).clone());
//            namePatches.add("");
            classPatches.add(-1);
        }
        inimg = inputImage;
        inimgarray = (byte[]) inimg.getProcessor().getPixels();
        pixelsMap = inputMap;
        D = inputD;
        steps=inpSteps;
        inSigma=inpSigma;
        fdsize=inpfdsize;
        fdbins=inpfdbins;
        min_size=inpmin_size;
        max_size=inpmax_size;
        upscale=inpupscale;
    }
public static void load(ImagePlus inputImage, ArrayList<int[]> inputLocXY, byte[] inputMap, int inputD) {
        siftFeatures.clear();
        locXY.clear();
        classPatches.clear();

        for (int i = 0; i < inputLocXY.size(); i++) {
            siftFeatures.add(new Vector<Feature>());
            locXY.add(inputLocXY.get(i).clone());
//            namePatches.add("");
            classPatches.add(-1);
        }
        inimg = inputImage;
        inimgarray = (byte[]) inimg.getProcessor().getPixels();
        pixelsMap = inputMap;
        D = inputD;
//        steps=inpSteps;
//        inSigma=inpSigma;
//        fdsize=inpfdsize;
//        fdbins=inpfdbins;
//        min_size=inpmin_size;
//        max_size=inpmax_size;
//        upscale=inpupscale;
    }

    public void run() {
        byte[] patchyeah = new byte[D * D];
        int W = inimg.getWidth();
        int x, y;
        for (int locIndex = begN; locIndex < endN; locIndex++) {
            //get a patch
            x = locXY.get(locIndex)[0];
            y = locXY.get(locIndex)[1];
            //obtain the sift features for the patch
            Vector<Feature> siftVector = calculateSIFT(inimgarray, W, x, y, D, patchyeah);
            siftFeatures.set(locIndex, siftVector);
            int pos = y * inimg.getWidth() + x;
            classPatches.set(locIndex, (pixelsMap[pos] == (byte) 255) ? 1 : 0);

        }
    }

//            rec.setRect(x, y, D, D);
//            Roi rec_roi = new Roi(rec);
//            ImageProcessor ipCopy = inimg.getChannelProcessor().duplicate();
//            ipCopy.setRoi(rec_roi);
//            ipCopy = ipCopy.crop();
//            ImagePlus impCopy = new ImagePlus("", ipCopy);
//            String namePatch = inimg.getShortTitle()
//                    + ",X,Y,D,i," + IJ.d2s(x, 0) + "," + IJ.d2s(y, 0) + "," + IJ.d2s(D, 0) + ".tif";
//            namePatches.set(locIndex, namePatch);
    //compare to know if it is a positive patch or not
//            int offset = y * inimg.getWidth();
//            int pos = offset + x;
//
//            if (pixelsMap[pos] == (byte) 255) //it's a neuron
//            {
//                classPatches.set(locIndex, 1); //add the class of the patch
////                                String filename = inimg.getShortTitle() + ",X,Y,D,i," + IJ.d2s(x, 0) + "," + IJ.d2s(y, 0) + "," + IJ.d2s(D, 0) + "," + IJ.d2s(count, 0) + ".tif";
////                                IJ.log(filename);
////                Color col = Color.YELLOW;
////                rec_roi.setFillColor(new Color(col.getRed(), col.getGreen(), col.getBlue(), 40));
////                ov.add(rec_roi);
////
////                inimg.setOverlay(ov);
////                inimg.setRoi(rec_roi);
////
////                inimg.updateAndDraw();
////                auxnn++;
////                nn++;
//            } else { //if it's noNeuron
//                classPatches.set(locIndex, 0);
//            }
    //to save the patches
//            if (save) {
//                //extract the patch
//                File f = new File(inputDirectory + "\\patches");
//                if (!f.exists()) {
//                    f.mkdirs();
//                }
//                fs = new FileSaver(impCopy);
//                String filename = f.getAbsolutePath() + File.separator + inimg.getShortTitle()
//                        + ",X,Y,D,i," + IJ.d2s(x, 0) + "," + IJ.d2s(y, 0) + "," + IJ.d2s(D, 0) + "," + IJ.d2s(count, 0) + ".tif";
//                fs.saveAsTiff(filename);
//            }
//            count++;
//                    if (false) {
    public static Vector<Feature> calculateSIFT(ImageProcessor imp) {
//        ImageProcessor ip1 = imp.getProcessor().convertToFloat();

        if (upscale) {
            scale = 2.0f;
        }
        Vector< Feature> fs1;
        FloatArray2DSIFT sift = new FloatArray2DSIFT(fdsize, fdbins);
        FloatArray2D fa = ImageArrayConverter.ImageToFloatArray2D(imp.convertToFloat());
        Filter.enhance(fa, 1.0f);

        if (upscale) {
            FloatArray2D fat = new FloatArray2D(fa.width * 2 - 1, fa.height * 2 - 1);
            FloatArray2DScaleOctave.upsample(fa, fat);
            fa = fat;
            fa = Filter.computeGaussianFastMirror(fa, (float) Math.sqrt(inSigma * inSigma - 1.0));
        } else {
            fa = Filter.computeGaussianFastMirror(fa, (float) Math.sqrt(inSigma * inSigma - 0.25));
        }
        sift.init(fa, steps, inSigma, min_size, max_size);
        fs1 = sift.run(max_size);
        Collections.sort(fs1);
        return fs1;
    }

    /*
     we'll try to replace this implementation with a new one without .duplicate for mosaics
     public static Vector<Feature> calculateSIFT(ImagePlus imp, int x, int y, int D) {
     Rectangle rec = new Rectangle();
     rec.setRect(x, y, D, D);
     Roi rec_roi = new Roi(rec);
     ImageProcessor ipCopy = imp.getChannelProcessor().duplicate(); // consumption!!!
     ipCopy.setRoi(rec_roi);
     ipCopy = ipCopy.crop();
     //        ImagePlus impCopy = new ImagePlus("", ipCopy);
     //        String namePatch = inimg.getShortTitle()
     //                + ",X,Y,D,i," + IJ.d2s(x, 0) + "," + IJ.d2s(y, 0) + "," + IJ.d2s(D, 0) + ".tif";
     //        namePatches.set(locIndex, namePatch);
     //obtain the sift features for the patch
     Vector<Feature> siftVector = calculateSIFT(ipCopy);
     return siftVector;
     }
     */
    public static Vector<Feature> calculateSIFT(byte[] inimg, int inimgW, int x, int y, int D, byte[] patch) {

        for (int k = 0; k < D * D; k++) {

            int i = k / D;
            int j = k % D;

            patch[(0 + j) * D + (0 + i)] = inimg[(y + j) * inimgW + (x + i)];

        }
        Vector<Feature> siftVector = calculateSIFT(new ByteProcessor(D, D, patch));
        return siftVector;

//            for (int i = 0; i < D; i++) {
//                for (int j = 0; j < D; j++) {
//                    patch[(0+j)*D+(0+i)] = inimg[(y+j)*inimgW+(x+i)];
//                }
//            }
//            ImagePlus t = new ImagePlus("", new ByteProcessor(D, D, patch));
//            imppatch.setProcessor();
    }

}
