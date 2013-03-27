package MLdetection;

import ij.*;
import ij.gui.*;
import ij.plugin.PlugIn;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.Image;
import imagescience.utility.I5DResource;
import java.awt.Button;
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
import java.util.ArrayList;
import javax.swing.JFileChooser;

/**
 *  Description of the Class
 *
 * @author     ismal
 * @created    June 2, 2006
 */
public class MLdetectionImage_ implements PlugIn, ActionListener, WindowListener {

private ImagePlus imp;
private int screenWidth,  screenHeight;
private String fileNameImage,  fileNamePositiveExamples,  fileNameTestExamples;
private FileWriter fw;
private Button changeImage,  changePositiveSamples,  changeTestSamples;
private GenericDialog dlg;
private int patchSize = 10;
private int adaT = 10;
private int weakStep;
private TextField fnImage,  fnTrainP,  fnTest;

public void run(String arg0) {

    // close all windows
    int[] c = WindowManager.getIDList();
    if (c != null) {
        for (int i = 0; i < c.length; i++) {
            WindowManager.getImage(c[i]).close();
        }
    }

    Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
    screenWidth = screenSize.width;
    screenHeight = screenSize.height;

    // selects the "point selection"-tool
    ij.IJ.setTool(7);

    // open config file with the filenames
    String plugins_name = IJ.getDirectory("plugins");
    String congifName = "mldetector.config";
    String confstr1 = " ", confstr2 = " ", confstr3 = " ";
    File file = new File(plugins_name + congifName);
    if (!file.exists()) {
        IJ.write("no config file mldetector.conf in plugins folder");
    } else {
        try {
            BufferedReader br = new BufferedReader(new FileReader(plugins_name + congifName));
            confstr1 = br.readLine();
            confstr2 = br.readLine();
            confstr3 = br.readLine();
            fileNameImage = br.readLine();
            fileNameTestExamples = br.readLine();
            fileNamePositiveExamples = br.readLine();
            patchSize = Integer.valueOf(br.readLine());
            adaT = Integer.valueOf(br.readLine());
            weakStep = Integer.valueOf(br.readLine());
            br.close();
        } catch (IOException e) {
            IJ.error("" + e);
            return;
        }
    }


    /* Build dialog: */
    dlg = new GenericDialog("Machine Learning - Detector");
    dlg.addPanel(makeButtonPanel_image(dlg));
    dlg.addPanel(makeButtonPanel_training_p(dlg));
    dlg.addPanel(makeButtonPanel_training_t(dlg));
    dlg.showDialog();


    if (dlg.wasCanceled()) {
        IJ.error("PlugIn canceled!");
        return;
    }

    fileNameImage = fnImage.getText();
    fileNamePositiveExamples = fnTrainP.getText();
    fileNameTestExamples = fnTest.getText();

    /* open image */
    IJ.open(fileNameImage);
    WindowManager.getCurrentWindow().addWindowListener(this);
    try {
        fw = new FileWriter(plugins_name + "mldetector.config");
        fw.write(confstr1 + "\n");
        fw.write(confstr2 + "\n");
        fw.write(confstr3 + "\n");
        fw.write(fileNameImage + "\n");
        fw.write(fileNameTestExamples + "\n");
        fw.write(fileNamePositiveExamples + "\n");
        fw.write(patchSize + "\n");
        fw.write(adaT + "\n");
        fw.write(weakStep + "\n");
        fw.close();
    } catch (IOException e) {
        IJ.error("Unable to save particle positions");
    }
    imp = ij.WindowManager.getCurrentImage();


    //################################## test the classifier ###################################

    //read adaboost coefficients
    double[][] adaboost;
    try {
        BufferedReader br = new BufferedReader(new FileReader(fileNamePositiveExamples));
        patchSize = Integer.valueOf(br.readLine());
        int leng = Integer.valueOf(br.readLine());
        adaboost = new double[leng][3];
        for (int i = 0; i < leng; i++) {
            adaboost[i][0] = Double.valueOf(br.readLine());
            adaboost[i][1] = Double.valueOf(br.readLine());
            adaboost[i][2] = Double.valueOf(br.readLine());
        }
        br.close();
    } catch (IOException e) {
        IJ.error("" + e);
        return;
    }

    // create the feature pool
    FeaturePool featurePool = new FeaturePool(patchSize);
    int fSize = featurePool.getSize();
    IJ.write("fearute pool size is " + fSize);
    //featurePool.visualizeFeaturePool();
    //ij.WindowManager.getCurrentWindow().setLocationAndSize(10, screenHeight - 250, 200, 200);


//    ArrayList<MMTrack> tracks_t = MdfFileParser.mdfFileRead(fileNameTestExamples);
//    int count = 0;
//    for (int i = 0; i < tracks_t.size(); i++) {
//        for (int j = 0; j < tracks_t.get(i).getSize(); j++) {
//            count++;
//        }
//    }
//
//    int nTest = count;
    Coordinates cin = new Coordinates();
    Coordinates cout = new Coordinates();
    Image inimg = Image.wrap(imp);

    inimg.axes(Axes.X + Axes.Y);
    Dimensions dims = inimg.dimensions();
    Image outimg = Image.create(dims, "imagescience.image.ByteImage");
    outimg.axes(Axes.X + Axes.Y);
//    Dimensions testdims = new Dimensions(patchSize, patchSize, 1, nTest, 1);
//    Image outimgT = Image.create(testdims, "imagescience.image.ByteImage");
//    outimgT.axes(Axes.X + Axes.Y);


//    String new_name = fileNameTestExamples.substring(0, fileNameTestExamples.length() - 4) + ".temp.mdf";
//    try {
//        fw = new FileWriter(new_name);
//        fw.write("MTrackJ 1.0.0 Data File\r\n");
//        fw.write("Displaying true true true true 0 1 0 100 18 1 0 0 0 1 true false true false\r\n");
//        fw.write("Assembly 1 FFFF00\r\n");
//        fw.write("Cluster 1 FFFF00\r\n");
//    } catch (IOException e) {
//        ij.IJ.error("Unable to save particle positions");
//    }

    //go through the whole image and check all the patches
    double[][] imageT = new double[patchSize][patchSize];
    for (cin.t = 0; cin.t < dims.t; cin.t++) {
        for (cin.y = 0; cin.y < dims.y - patchSize; cin.y++) {
            for (cin.x = 0; cin.x < dims.x - patchSize; cin.x++) {
                inimg.get(cin, imageT);
                int[][] imageIntT = getIntegralImage(imageT);
                int[] imFeaturesT = computeSelectedFeatures(imageIntT, featurePool, adaboost);
                int test = applyAdaBoost(adaboost, imFeaturesT);
                if (test == 1) {
                    cout.x = cin.x + patchSize / 2;
                    cout.y = cin.y + patchSize / 2;
                    cout.t = cin.t;
                    outimg.set(cout, 255);
                } else {
                }
            }
        }
    }
//
//    count = 0;
//    int count_p = 0;
//    int count_n = 0;
//    String result = "";
//    for (int i = 0; i < nTest; i++) {
//        try {
//            fw.write("Track " + (i + 1) + " 00FF00\r\n");
//        } catch (IOException e) {
//            IJ.error("" + e);
//        }
//        for (int j = 0; j < tracks_t.get(i).getSize(); j++) {
//            cin.x = Math.round(tracks_t.get(i).getX(j) - patchSize / 2);
//            cin.y = Math.round(tracks_t.get(i).getY(j) - patchSize / 2);
//            cin.t = Math.round(tracks_t.get(i).gett(j) - 1);
//            inimg.get(cin, imageT);
//
//
//            //double[][] interpolatedImage = getInterpolatedImage(imageT);
//            //int[][] imageIntT = getIntegralImage(interpolatedImage);
//            int[][] imageIntT = getIntegralImage(imageT);
//
//            int[] imFeaturesT = computeSelectedFeatures(imageIntT, featurePool, adaboost);
//            int test = applyAdaBoost(adaboost, imFeaturesT);
//            if (test == 1) {
//                count_p++;
//            } else {
//                count_n++;
//            }
//
//            try {
//                fw.write("Point " + (j + 1) +
//                         " " + IJ.d2s(tracks_t.get(i).getX(j), 2) +
//                         " " + IJ.d2s(tracks_t.get(i).getY(j), 2) +
//                         " " + IJ.d2s(tracks_t.get(i).getZ(j), 2) +
//                         " " + IJ.d2s(tracks_t.get(i).gett(j), 2) +
//                         " " + IJ.d2s(2 - test, 2) + "\r\n");
//            } catch (IOException e) {
//                IJ.error("" + e);
//            }
//
//            result += test + " ";
//
//            cout.t = count;
//            outimgT.set(cout, imageT);
//
//            count++;
//        }
//    }
//
//    try {
//        fw.write("End of MTrackJ Data File\r\n");
//        fw.close();
//
//
//    } catch (IOException e) {
//        IJ.error("Unable to save particle positions");
//    }
//    IJ.write("positive responses : " + count_p + " , negative : " + count_n + "   " + (count_p / 4096.0) + "   " + (count_n / 4096.0));
//
//    Image testimage = Image.create(new Dimensions(imp.getWidth(), imp.getHeight(), 1, imp.getNFrames(), 2),
//                                   "imagescience.image.ByteImage");
//    Coordinates ctest = new Coordinates(0, 0, 0, 0, 0);
//    testimage.axes(Axes.X + Axes.Y);
//    double[][] test = new double[imp.getHeight()][imp.getWidth()];
//
//    Image test1 = Image.wrap(imp);
//    test1.axes(Axes.X + Axes.Y);
//    for (int i = 0; i < imp.getNFrames(); i++) {
//        ctest.t = i;
//        test1.get(ctest, test);
//        ctest.s = 0;
//        testimage.set(ctest, test);
//        ctest.s = 1;
//        testimage.set(ctest, test);
//        ctest.s = 0;
//    }
        ImagePlus outimp = outimg.imageplus();
        boolean i5dinstalled = false;
        try {
            Class.forName("i5d.Image5D");
            i5dinstalled = true;
        } catch (Exception e) {
        }
        if (i5dinstalled) {
            final Dimensions outdims = outimg.dimensions();
            if (outdims.t > 1 || outdims.c > 1) {
                outimp = I5DResource.convert(outimp, true);
                //outimp.setDisplayMode(2);
                if (I5DResource.instance(imp)) {
                    I5DResource.transfer(imp, outimp);
                }
            }
        }
        //ImagePlus outimp = outimg.imageplus();
        outimp.setTitle("Output");
        outimp.show();

//    IJ.runPlugIn("MTrackJ_", new_name);
        //IJ.write(result);
//        ImagePlus outimpT = outimgT.imageplus();
//        outimpT.show();
//        for (int i = 0; i < 7; i++) {
//            outimpT.getCanvas().zoomIn(0, 0);
//        }
    }

private

 int applyAdaBoost(double[][] adaboost, int[] imFeaturesT) {
    double res = 0;
    double object_res = 0;
    for (int i = 0; i < adaboost.length; i++) {
        object_res += adaboost[i][1];
    }
    object_res *= 0.5;

    for (int i = 0; i < adaboost.length; i++) {
        if (applyClassifier(imFeaturesT[i], adaboost[i][2])) {
            res += adaboost[i][1];
        }
    }
    int test = (res > object_res) ? 1 : 0;

    return test;
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

private int[] computeSelectedFeatures(int[][] imageInt, FeaturePool featurePool, double[][] adaboost) {
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

/** 
 * Creates a panel containing change image and change neurofile buttons
 */
Panel makeButtonPanel_image(GenericDialog gd) {
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

    changePositiveSamples = new Button("AdaBoost coefs.");
    changePositiveSamples.addActionListener(this);
    buttons.add(changePositiveSamples);

    fnTrainP = new TextField(fileNamePositiveExamples, 80);
    fnTrainP.addActionListener(this);
    buttons.add(fnTrainP);

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