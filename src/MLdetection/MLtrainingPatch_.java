package MLdetection;

import ij.*;
import ij.gui.*;
import ij.plugin.PlugIn;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.Image;
import imagescience.utility.Progressor;
import java.awt.Button;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Panel;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import javax.swing.JFileChooser;

/**
 *  Description of the Class
 *
 * @author     ismal
 * @created    June 2, 2006
 */
public class MLtrainingPatch_ implements PlugIn, ActionListener, WindowListener {

private ImagePlus imp;
private int screenWidth,  screenHeight;
private String fileNameImage,  fileNamePositiveExamples,  fileNameNegativeExamples,  fileNameTestExamples;
private FileWriter fw;
private Button changeImage,  changePositiveSamples,  changeNegativeSamples,  changeTestSamples;
private GenericDialog dlg;
private int patchSize = 10;
private int weakStep = 1000;
private int adaT = 10;
private TextField fnImage,  fnTrainP,  fnTrainN,  fnTest;
Progressor pgs;

public void run(String arg0) {

    // close all windows
    int[] c = WindowManager.getIDList();
    if (c != null) {
        for (int i = 0; i < c.length; i++) {
            WindowManager.getImage(c[i]).close();
        }
    }

    // get the desktop size in order to place the "Result" and image window correctly
    Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
    screenWidth = screenSize.width;
    screenHeight = screenSize.height;

    // selects the "point selection"-tool
    ij.IJ.setTool(7);

    ////////////////////////////////////////////////////////////////////////////
    // open config file with the filenames
    // the config file is
    // line 1 - image for training
    // line 2 - positiove samples for that image
    // line 3 - negative samples for that image 
    // line 4 - last test/classification image 
    // line 5 - .mdf file with point for the classification 
    // line 6 = .adaboost file with coefficients
    // line 7 - Size of the scanning window
    // line 8 - Number of adaboost runs 
    // line 9 = Search Step for weak classifier 

    // !!!! this plugin modifies only the first three lines of the config file, 
    // !!!! the rest of the file is simply rewritten

    String plugins_name = IJ.getDirectory("plugins");
    String congifName = "mldetector.config";
    // dummy string for rewriting the unrelated strings from the config file
    String confstr4 = " ", confstr5 = " ", confstr6 = " ";
    File file = new File(plugins_name + congifName);
    if (!file.exists()) {
        fileNameImage = "...select the image file...";
        fileNamePositiveExamples = "...select the mdf file with positive samples...";
        fileNameNegativeExamples = "...select the mdf file with negative samples...";
    } else {
        try {
            BufferedReader br = new BufferedReader(new FileReader(plugins_name + congifName));
            fileNameImage = br.readLine();
            fileNamePositiveExamples = br.readLine();
            fileNameNegativeExamples = br.readLine();
            confstr4 = br.readLine();
            confstr5 = br.readLine();
            confstr6 = br.readLine();
            patchSize = Integer.valueOf(br.readLine());
            adaT = Integer.valueOf(br.readLine());
            weakStep = Integer.valueOf(br.readLine());
            br.close();
        } catch (IOException e) {
            IJ.error("" + e);
            return;
        }
    }

    Prefs.get("miro.adaNumber", 2.0);
    Prefs.get("miro.path", IJ.getDirectory("plugins"));
     
    /* Build dialog: */
    dlg = new GenericDialog("Machine Learning : Training");
    //dlg.addNumericField("Size of the scanning window : ", patchSize, 0);
    dlg.addNumericField("Number of adaboost runs : ", adaT, 0);
    dlg.addNumericField("Search Step for weak classifier : ", weakStep, 0);
    //dlg.addPanel(makeButtonPanel_image(dlg));
    dlg.addPanel(makeButtonPanel_training_p(dlg));
    dlg.addPanel(makeButtonPanel_training_n(dlg));
    dlg.showDialog();
    if (dlg.wasCanceled()) {
        IJ.error("PlugIn canceled!");
        return;
    }

   // patchSize = (int) dlg.getNextNumber();
    adaT = (int) dlg.getNextNumber();
    weakStep = (int) dlg.getNextNumber();

    Prefs.set("miro.adaNumber", adaT);
    
    //fileNameImage = fnImage.getText();
    fileNamePositiveExamples = fnTrainP.getText();
    fileNameNegativeExamples = fnTrainN.getText();

    try {
        fw = new FileWriter(plugins_name + "mldetector.config");
        fw.write(fileNameImage + "\n");
        fw.write(fileNamePositiveExamples + "\n");
        fw.write(fileNameNegativeExamples + "\n");
        fw.write(confstr4 + "\n");
        fw.write(confstr5 + "\n");
        fw.write(confstr6 + "\n");
        fw.write(patchSize + "\n");
        fw.write(adaT + "\n");
        fw.write(weakStep + "\n");
        fw.close();
    } catch (IOException e) {
        IJ.error("Unable to save particle positions");
    }
    
    /* open image */
    IJ.open(fileNamePositiveExamples);
    WindowManager.getCurrentWindow().addWindowListener(this);
    imp = ij.WindowManager.getCurrentImage();
    patchSize = imp.getWidth();
    for (int i = 0; i < 7; i++) {
        imp.getCanvas().zoomIn(0, 0);
    }
    imp.getWindow().setLocationAndSize(160, screenHeight - 350, 200, 200);

    // create the feature pool 
    FeaturePool featurePool = new FeaturePool(patchSize);
    int fSize = featurePool.getSize();
    IJ.log("the fearute pool of size " + fSize + " is created");
    featurePool.visualizeFeaturePool();
    ij.WindowManager.getCurrentWindow().setLocationAndSize(10, screenHeight - 350, 200, 200);

    ////////////////////////////////////////////////////////////////////////////
    // extract the positive samples
    ////////////////////////////////////////////////////////////////////////////
    Coordinates cin = new Coordinates();
    Image inimg = Image.wrap(imp);
    inimg.axes(Axes.X + Axes.Y);
    Dimensions testdims = inimg.dimensions();

    double[][] imageP = new double[patchSize][patchSize];
    int[][] imFeaturesP = new int[testdims.t][];
    for (int i = 0; i < testdims.t; i++) {
            cin.t = i;
            inimg.get(cin, imageP);
            int[][] imageIntP = getIntegralImage(imageP);
            imFeaturesP[i] = computeFeatures(imageIntP, featurePool);
    }

    ////////////////////////////////////////////////////////////////////////////
    // extract the negative samples
    ////////////////////////////////////////////////////////////////////////////

    /* open image */
    IJ.open(fileNameNegativeExamples);
    WindowManager.getCurrentWindow().addWindowListener(this);
    imp = ij.WindowManager.getCurrentImage();
    patchSize = imp.getWidth();
    for (int i = 0; i < 7; i++) {
        imp.getCanvas().zoomIn(0, 0);
    }
    imp.getWindow().setLocationAndSize(310, screenHeight - 350, 200, 200);

    inimg = Image.wrap(imp);
    inimg.axes(Axes.X + Axes.Y);
    testdims = inimg.dimensions();
    
    double[][] imageN = new double[patchSize][patchSize];
    int[][] imFeaturesN = new int[testdims.t][];
    for (int i = 0; i < testdims.t; i++) {
            cin.t = i;
            inimg.get(cin, imageN);
            int[][] imageIntN = getIntegralImage(imageN);
            imFeaturesN[i] = computeFeatures(imageIntN, featurePool);
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // train the classifier
    ////////////////////////////////////////////////////////////////////////////

    double[][] adaboost = trueAdaBoost(imFeaturesP, imFeaturesN, adaT);

    // save the coefio file 
    try {
        fw = new FileWriter(fileNamePositiveExamples + "." + adaT + ".adaboost");
        fw.write(IJ.d2s(patchSize, 0) + "\r\n");
        fw.write(IJ.d2s(adaboost.length, 0) + "\r\n");
        for (int i = 0; i < adaboost.length; i++) {
            fw.write(IJ.d2s(adaboost[i][0], 0) + "\r\n");
            fw.write(IJ.d2s(adaboost[i][1], 6) + "\r\n");
            fw.write(IJ.d2s(adaboost[i][2], 6) + "\r\n");
        }
        fw.close();
    } catch (IOException e) {
        ij.IJ.error("Unable to save particle positions");
    }


    // plot the features with thresholds 
    ImagePlus imp1 = new ImagePlus();
    ImageStack imstack = new ImageStack(528, 255);
    for (int i = 0; i < adaboost.length; i++) {
        Plot plot = plotFearutePerAllSamples(imFeaturesP, imFeaturesN, (int) adaboost[i][0], adaboost[i][2]);
        //Plot plot = plotFearutePerAllSamples(imFeaturesP, imFeaturesN, i, 0);
        imstack.addSlice("thresh = " + IJ.d2s(adaboost[i][2], 3), plot.getProcessor());
    }
    imp1.setStack("test", imstack);
    imp1.show();

}

private double[][] trueAdaBoost(int[][] imFeaturesP, int[][] imFeaturesN, int T) {
    int sizeP = imFeaturesP.length;
    int sizeN = imFeaturesN.length;
    int fSize = imFeaturesN[0].length;
    int m = sizeP + sizeN;
    double[] w = new double[m];
    double[][] adaboost = new double[T][3]; // [id, alpha, threshold]
    //initial weights

    double q = 1.0 / m;
    for (int i = 0; i < m; i++) {
        w[i] = q;
    }
    int tt = T;
    for (int t = 0; t < T; t++) {
        double[][] thresh = weightedWeakClassification(imFeaturesP, imFeaturesN, w, weakStep); // [feature_i][optimal threshold and the score /error]
        //find minimum error 
        double mineps = thresh[0][1];
        int bestClassifier = 0;
        for (int k = 1; k < fSize; k++) {
            if (mineps > thresh[k][1]) {
                mineps = thresh[k][1];
                bestClassifier = k;
            }
        }
        adaboost[t][0] = bestClassifier;
        if (mineps == 0 || mineps > 0.5) {
            adaboost[t][1] = 1;
            adaboost[t][0] = bestClassifier;
            adaboost[t][2] = thresh[bestClassifier][0];

            IJ.log(t + " min eps = " + mineps + "   " +
                    bestClassifier + "   " + thresh[bestClassifier][0]);
            tt = t;
            break;
        }
        double beta = mineps / (1 - mineps);

        //update weights
        for (int i = 0; i < sizeP; i++) {
            if (applyClassifier(imFeaturesP[i][bestClassifier], thresh[bestClassifier][0])) {
                w[i] *= beta;
            }
        }
        for (int i = 0; i < sizeN; i++) {
            if (!applyClassifier(imFeaturesN[i][bestClassifier], thresh[bestClassifier][0])) {
                w[i + sizeP] *= beta;
            }
        }
        adaboost[t][1] = Math.log(1 / beta);
        adaboost[t][2] = thresh[bestClassifier][0];

        //normalize weights
        double sum = 0;
        for (int i = 0; i < m; i++) {
            sum += w[i];
        }
        for (int i = 0; i < m; i++) {
            w[i] /= sum;
        }
        // number of runs - id of the feature - alpha - threshold - error 
        IJ.log("t = " + (t + 1) + " -> best classifier:" + IJ.d2s(adaboost[t][0] + 1, 0) + " Math.log(1 / beta):  " +
                IJ.d2s(adaboost[t][1], 4) + "  optimal threshold: " + IJ.d2s(adaboost[t][2], 4) + "  error:   " +
                IJ.d2s(mineps, 6));
    }
    if (tt != T) {
        double[][] new_adaboost = new double[tt + 1][3];
        for (int i = 0; i < tt + 1; i++) {
            new_adaboost[i] = adaboost[i];
        }
        return new_adaboost;
    } else {
        return adaboost;
    }
}

private double[][] weightedWeakClassification(int[][] imFeaturesP, int[][] imFeaturesN,
        double[] w, int dt) {
    int sizeP = imFeaturesP.length;
    int sizeN = imFeaturesN.length;
    int fSize = imFeaturesN[0].length;
    double[][] thresh = new double[fSize][2];

    for (int k = 0; k < fSize; k++) {
        // find max and min
        double max = imFeaturesP[0][k];
        double min = imFeaturesP[0][k];
        for (int i = 1; i < sizeP; i++) {
            max = (max < imFeaturesP[i][k]) ? imFeaturesP[i][k] : max;
            min = (min > imFeaturesP[i][k]) ? imFeaturesP[i][k] : min;
        }
        for (int i = 0; i < sizeN; i++) {
            max = (max < imFeaturesN[i][k]) ? imFeaturesN[i][k] : max;
            min = (min > imFeaturesN[i][k]) ? imFeaturesN[i][k] : min;
        }

        double step = (max - min) / dt;
        double thr = min;
        double[][] count = new double[dt][2];
        for (int j = 0; j < dt - 1; j++) {
            thr += step;
            double score = 0.0;
            for (int i = 0; i < sizeP; i++) {
                if (!applyClassifier(imFeaturesP[i][k], thr)) {
                    score += w[i];
                }
            }
            for (int i = 0; i < sizeN; i++) {
                if (applyClassifier(imFeaturesN[i][k], thr)) {
                    score += w[i + sizeP];
                }
            }
            count[j][1] = score;
            count[j][0] = thr;
        }
        // find optimal threshold
        thresh[k][1] = count[0][1];
        for (int j = 1; j < dt - 1; j++) {
            if (thresh[k][1] > count[j][1]) {
                thresh[k][1] = count[j][1];
                thresh[k][0] = count[j][0];
            }
        }
    }
//    // show all features
//    ImagePlus imp1 = new ImagePlus();
//    ImageStack imstack = new ImageStack(528, 255);
//    //imp1.setDimensions(1, 1, adaboost.length);
//    for (int i = 0; i < imFeaturesN[0].length; i++) {
//        //Plot plot = plotFearutePerAllSamples(imFeaturesP, imFeaturesN, (int)adaboost[i][0], adaboost[i][2]);
//        Plot plot = plotFearutePerAllSamples(imFeaturesP, imFeaturesN, i, thresh[i][0]);
//        imstack.addSlice("thresh = " + thresh[i][0], plot.getProcessor());
//    }
//    imp1.setStack("test", imstack);
//    imp1.show();

    return thresh;
}

private boolean applyClassifier(double x, double thresh) {
    return (x >= thresh) ? true : false;
}

private int[][] getIntegralImage(double[][] imageZ) {
    int dimsx = imageZ[0].length;
    int dimsy = imageZ.length;
    int[][] imageInt = new int[dimsy][dimsx];
    int[][] s = new int[dimsy][dimsx];
    for (int j = 0; j < dimsy; j++) {
        for (int i = 0; i < dimsx; i++) {
            s[j][i] = (j - 1 >= 0) ? s[j - 1][i] + (int) imageZ[j][i] : (int) imageZ[j][i];
            imageInt[j][i] = (i - 1 >= 0) ? imageInt[j][i - 1] + s[j][i] : s[j][i];
        }
    }
    return imageInt;
}

private double[][] getInterpolatedImage(double[][] imageZ) {
    int dimsx = imageZ[0].length;
    int dimsy = imageZ.length;
    double[][] imageInt = new double[2 * dimsy][2 * dimsx];
    for (int j = 0; j < dimsy; j++) {
        for (int i = 0; i < dimsx; i++) {
            imageInt[2 * j][2 * i] = imageZ[j][i];
            imageInt[2 * j + 1][2 * i] = imageZ[j][i];
            imageInt[2 * j][2 * i + 1] = imageZ[j][i];
            imageInt[2 * j + 1][2 * i + 1] = imageZ[j][i];
        }
    }
    return imageInt;
}

private int[] computeFeatures(int[][] imageInt, FeaturePool featurePool) {
    int dimsx = imageInt[0].length;
    int dimsy = imageInt.length;
    int fSize = featurePool.getSize();
    int[] imFeatures = new int[fSize];
    for (int i = 0; i < fSize; i++) {
        int[] r0 = featurePool.getFeature(i).getWhiteSquare();
        int[] r1 = featurePool.getFeature(i).getBlackSquare();
        double[] w = featurePool.getFeature(i).getWeights();

        int i11 = imageInt[r1[1] + r1[3] - 1][r1[0] + r1[2] - 1];
        int i00 = ((r1[0] - 1 < 0) || (r1[1] - 1 < 0)) ? 0 : imageInt[r1[1] - 1][r1[0] - 1];
        int i01 = ((r1[0] - 1 < 0)) ? 0 : imageInt[r1[1] + r1[3] - 1][r1[0] - 1];
        int i10 = ((r1[1] - 1 < 0)) ? 0 : imageInt[r1[1] - 1][r1[0] + r1[2] - 1];
        int sum1 = i11 + i00 - i10 - i01;

        i11 = imageInt[r0[1] + r0[3] - 1][r0[0] + r0[2] - 1];
        i00 = ((r0[0] - 1 < 0) || (r0[1] - 1 < 0)) ? 0 : imageInt[r0[1] - 1][r0[0] - 1];
        i01 = ((r0[0] - 1 < 0)) ? 0 : imageInt[r0[1] + r0[3] - 1][r0[0] - 1];
        i10 = ((r0[1] - 1 < 0)) ? 0 : imageInt[r0[1] - 1][r0[0] + r0[2] - 1];
        int sum0 = i11 + i00 - i10 - i01;

        imFeatures[i] = (int) (w[0] * sum0 + w[1] * sum1);
    }

    return imFeatures;
}

private int[] computeSelectedFeatures(int[][] imageInt, FeaturePool featurePool,
        double[][] adaboost) {
    int dimsx = imageInt[0].length;
    int dimsy = imageInt.length;
    int size = adaboost.length;
    int[] imFeatures = new int[size];
    for (int i = 0; i < size; i++) {
        int j = (int) adaboost[i][0];
        int[] r0 = featurePool.getFeature(j).getWhiteSquare();
        int[] r1 = featurePool.getFeature(j).getBlackSquare();
        double[] w = featurePool.getFeature(j).getWeights();

        int i11 = imageInt[r1[1] + r1[3] - 1][r1[0] + r1[2] - 1];
        int i00 = ((r1[0] - 1 < 0) || (r1[1] - 1 < 0)) ? 0 : imageInt[r1[1] - 1][r1[0] - 1];
        int i01 = ((r1[0] - 1 < 0)) ? 0 : imageInt[r1[1] + r1[3] - 1][r1[0] - 1];
        int i10 = ((r1[1] - 1 < 0)) ? 0 : imageInt[r1[1] - 1][r1[0] + r1[2] - 1];
        int sum1 = i11 + i00 - i10 - i01;

        i11 = imageInt[r0[1] + r0[3] - 1][r0[0] + r0[2] - 1];
        i00 = ((r0[0] - 1 < 0) || (r0[1] - 1 < 0)) ? 0 : imageInt[r0[1] - 1][r0[0] - 1];
        i01 = ((r0[0] - 1 < 0)) ? 0 : imageInt[r0[1] + r0[3] - 1][r0[0] - 1];
        i10 = ((r0[1] - 1 < 0)) ? 0 : imageInt[r0[1] - 1][r0[0] + r0[2] - 1];
        int sum0 = i11 + i00 - i10 - i01;

        imFeatures[i] = (int) (w[0] * sum0 + w[1] * sum1);
    }

    return imFeatures;
}

private Plot plotFearutePerAllSamples(int[][] imFeaturesP, int[][] imFeaturesN, int c, double thresh) {
    int size = imFeaturesP[0].length;
    int sizep = imFeaturesP.length;
    int sizen = imFeaturesN.length;
    double[] xp = new double[sizep];
    double[] xn = new double[sizen];
    double[] yp = new double[sizep];
    double[] yn = new double[sizen];
    int maxyp = imFeaturesP[0][c];
    int minyp = imFeaturesP[0][c];
    for (int i = 0; i < sizep; i++) {
        xp[i] = i;
        yp[i] = imFeaturesP[i][c];
        maxyp = (maxyp < imFeaturesP[i][c]) ? imFeaturesP[i][c] : maxyp;
        minyp = (minyp > imFeaturesP[i][c]) ? imFeaturesP[i][c] : minyp;
    }
    int maxyn = imFeaturesN[0][c];
    int minyn = imFeaturesN[0][c];
    for (int i = 0; i < sizen; i++) {
        xn[i] = i;
        yn[i] = imFeaturesN[i][c];
        maxyn = (maxyn < imFeaturesN[i][c]) ? imFeaturesN[i][c] : maxyn;
        minyn = (minyn > imFeaturesN[i][c]) ? imFeaturesN[i][c] : minyn;
    }

    Plot plot = new Plot("Average velocity = " + c, "Sample id", "Feature", xp, yp);
    plot.setLimits(0, Math.max(sizen, sizep), Math.min(minyn, minyp), Math.max(maxyn, maxyp));
    plot.setColor(Color.RED);
    plot.addPoints(xn, yn, Plot.LINE);
    plot.setColor(Color.BLACK);
    double[] linex = new double[]{0, Math.max(sizen, sizep)};
    double[] liney = new double[]{thresh, thresh};
    plot.setLineWidth(5);
    plot.addPoints(linex, liney, Plot.LINE);
    plot.setColor(Color.BLUE);
    plot.setLineWidth(2);
    //plot.show();
    return plot;
}


//////////////////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////// 
//  Creates a panel containing change image and change neurofile buttons
//////////////////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////// 

Panel makeButtonPanel_image() { //GenericDialog gd
    Panel buttons = new Panel();
    buttons.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 1));

    changeImage = new Button("Image....");
    changeImage.addActionListener(this);
    buttons.add(changeImage);

    fnImage = new TextField(fileNameImage, 80);
    fnImage.addActionListener(this);
    buttons.add(fnImage);

    return buttons;
}

Panel makeButtonPanel_training_p(GenericDialog gd) {
    Panel buttons = new Panel();
    buttons.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 1));

    changePositiveSamples = new Button("Positive set");
    changePositiveSamples.addActionListener(this);
    buttons.add(changePositiveSamples);

    fnTrainP = new TextField(fileNamePositiveExamples, 80);
    fnTrainP.addActionListener(this);
    buttons.add(fnTrainP);

    return buttons;
}

Panel makeButtonPanel_training_n(GenericDialog gd) {
    Panel buttons = new Panel();
    buttons.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 1));

    changeNegativeSamples = new Button("Negative set");
    changeNegativeSamples.addActionListener(this);
    buttons.add(changeNegativeSamples);

    fnTrainN = new TextField(fileNameNegativeExamples, 80);
    fnTrainN.addActionListener(this);
    buttons.add(fnTrainN);

    return buttons;
}

Panel makeButtonPanel_training_t(GenericDialog gd) {
    Panel buttons = new Panel();
    buttons.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 1));

    changeTestSamples = new Button("Test set");
    changeTestSamples.addActionListener(this);
    buttons.add(changeTestSamples);

    fnTest = new TextField(fileNameTestExamples, 80);
    fnTest.addActionListener(this);
    buttons.add(fnTest);

    return buttons;
}

public void actionPerformed(ActionEvent e) {
    Object source = e.getSource();
    if (source == changeImage) {
        JFileChooser fc = new JFileChooser(fileNameImage);
        //fc.addChoosableFileFilter(new MyFilterLSM());
        int returnVal = fc.showOpenDialog(null);
        if (returnVal != JFileChooser.APPROVE_OPTION) {
            return;
        }
        String fileName = fc.getSelectedFile().getAbsolutePath();
        IJ.showStatus("Opening: " + fileName);
        fnImage.setText(fileName);
    } else if (source == changePositiveSamples) {
        JFileChooser fc = new JFileChooser(fileNamePositiveExamples);
        //fc.addChoosableFileFilter(new MyFilterTXT());
        int returnVal = fc.showOpenDialog(null);
        if (returnVal != JFileChooser.APPROVE_OPTION) {
            return;
        }
        String fileName = fc.getSelectedFile().getAbsolutePath();
        IJ.showStatus("Opening: " + fileName);
        fnTrainP.setText(fileName);
    } else if (source == changeNegativeSamples) {
        JFileChooser fc = new JFileChooser(fileNameNegativeExamples);
        //fc.addChoosableFileFilter(new MyFilterTXT());
        int returnVal = fc.showOpenDialog(null);
        if (returnVal != JFileChooser.APPROVE_OPTION) {
            return;
        }
        String fileName = fc.getSelectedFile().getAbsolutePath();
        IJ.showStatus("Opening: " + fileName);
        fnTrainN.setText(fileName);

    } else if (source == changeTestSamples) {
        JFileChooser fc = new JFileChooser(fileNameTestExamples);
        //fc.addChoosableFileFilter(new MyFilterTXT());
        int returnVal = fc.showOpenDialog(null);
        if (returnVal != JFileChooser.APPROVE_OPTION) {
            return;
        }
        String fileName = fc.getSelectedFile().getAbsolutePath();
        IJ.showStatus("Opening: " + fileName);
        fnTest.setText(fileName);
    }
    dlg.repaint();
}

public void windowOpened(WindowEvent e) {
}

public void windowClosing(WindowEvent e) {
    int[] c = WindowManager.getIDList();
    for (int i = 0; i < c.length; i++) {
        WindowManager.getImage(c[i]).close();
    }
}

public void windowClosed(WindowEvent e) {
}

public void windowIconified(WindowEvent e) {
}

public void windowDeiconified(WindowEvent e) {
}

public void windowActivated(WindowEvent e) {
}

public void windowDeactivated(WindowEvent e) {
}
}