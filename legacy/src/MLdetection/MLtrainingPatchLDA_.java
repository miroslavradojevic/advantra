package MLdetection;

import Jama.Matrix;
import ij.*;
import ij.gui.*;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.Image;
import imagescience.utility.Progressor;
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
import javax.swing.JFileChooser;

/**
 *  Description of the Class
 *
 * @author     ismal
 * @created    June 2, 2006
 */
public class MLtrainingPatchLDA_ implements PlugIn, ActionListener, WindowListener {

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
        fileNamePositiveExamples = "...select the mdf fiel with positive samples...";
        fileNameNegativeExamples = "...select the mdf fiel with negative samples...";
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

    /* Build dialog: */
    dlg = new GenericDialog("Machine Leatning : Training LDA");
    //dlg.addNumericField("Size of the scanning window : ", patchSize, 0);
    //dlg.addNumericField("Number of adaboost runs : ", adaT, 0);
    //dlg.addNumericField("Search Step for weak classifier : ", weakStep, 0);
    //dlg.addPanel(makeButtonPanel_image(dlg));
    dlg.addPanel(makeButtonPanel_training_p(dlg));
    dlg.addPanel(makeButtonPanel_training_n(dlg));
    dlg.showDialog();
    if (dlg.wasCanceled()) {
        IJ.error("PlugIn canceled!");
        return;
    }

    //patchSize = (int) dlg.getNextNumber();
    //adaT = (int) dlg.getNextNumber();
    //weakStep = (int) dlg.getNextNumber();

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

    ////////////////////////////////////////////////////////////////////////////
    // extract the positive samples
    ////////////////////////////////////////////////////////////////////////////
    Coordinates cin = new Coordinates();
    Image inimg = Image.wrap(imp);
    inimg.axes(Axes.X + Axes.Y);
    Dimensions testdims = inimg.dimensions();
    double[][] arrPositive = new double[testdims.t][patchSize * patchSize];
    double[][] imageP = new double[patchSize][patchSize];
    for (int i = 0; i < testdims.t; i++) {
        cin.t = i;
        inimg.get(cin, imageP);
        for (int jj = 0; jj < patchSize; jj++) {
            for (int ii = 0; ii < patchSize; ii++) {
                arrPositive[i][jj * patchSize + ii] = imageP[jj][ii];
            }
        }
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
    double[][] arrNegative = new double[testdims.t][patchSize * patchSize];
    double[][] imageN = new double[patchSize][patchSize];
    for (int i = 0; i < testdims.t; i++) {
        cin.t = i;
        inimg.get(cin, imageN);
        for (int jj = 0; jj < patchSize; jj++) {
            for (int ii = 0; ii < patchSize; ii++) {
                arrNegative[i][jj * patchSize + ii] = imageN[jj][ii];
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // train the clssifier
    ////////////////////////////////////////////////////////////////////////////

    //compute means 
    int patchSize2 = patchSize * patchSize;
    Matrix mat_mean_fore = new Matrix(patchSize2, 1);
    Matrix mat_mean_back = new Matrix(patchSize2, 1);
    double[] mean_fore = new double[patchSize2];
    double[] mean_back = new double[patchSize2];
    for (int i = 0; i < patchSize2; i++) {
        for (int j = 0; j < arrPositive.length; j++) {
            mean_fore[i] += arrPositive[j][i];
        }
        mean_fore[i] /= arrPositive.length;
        mat_mean_fore.set(i, 0, mean_fore[i]);

        for (int j = 0; j < arrNegative.length; j++) {
            mean_back[i] += arrNegative[j][i];
        }
        mean_back[i] /= arrNegative.length;
        mat_mean_back.set(i, 0, mean_back[i]);
    }

    //compute variance 
    Matrix Sigma_fore = new Matrix(patchSize2, patchSize2);
    Matrix Sigma_back = new Matrix(patchSize2, patchSize2);
    for (int i = 0; i < patchSize2; i++) {
        for (int j = 0; j < patchSize2; j++) {
            double sum = 0;
            for (int k = 0; k < arrPositive.length; k++) {
                sum += (arrPositive[k][i] - mean_fore[i]) * (arrPositive[k][j] - mean_fore[j]);
            }
            sum /= arrPositive.length;
            Sigma_fore.set(i, j, sum);
            Sigma_fore.set(j, i, sum);

            sum = 0;
            for (int k = 0; k < arrNegative.length; k++) {
                sum += (arrNegative[k][i] - mean_back[i]) * (arrNegative[k][j] - mean_back[j]);
            }
            sum /= arrNegative.length;
            Sigma_back.set(i, j, sum);
            Sigma_back.set(j, i, sum);
        }
    }

    Matrix w = ((Sigma_back.plus(Sigma_fore)).inverse()).times(mat_mean_fore.minus(mat_mean_back));

    Matrix m_fore = w.transpose().times(mat_mean_fore);
    Matrix m_back = w.transpose().times(mat_mean_back);

    // save the coefio file 
    try {
        fw = new FileWriter(fileNamePositiveExamples + ".lda");
        fw.write(IJ.d2s(patchSize, 0) + "\r\n");
        for (int i = 0; i < patchSize2; i++) {
            fw.write(IJ.d2s(w.get(i, 0), 6) + "\r\n");
        }
        fw.write(IJ.d2s(m_fore.get(0, 0), 6) + "\r\n");
        fw.write(IJ.d2s(m_back.get(0, 0), 6) + "\r\n");
        fw.close();
    } catch (IOException e) {
        ij.IJ.error("Unable to save particle positions");
    }


    // plot the features with thresholds 
    ImagePlus imp1 = new ImagePlus();
    float[][] imageW = new float[patchSize][patchSize];
    for (int j = 0; j < patchSize; j++) {
        for (int i = 0; i < patchSize; i++) {
            imageW[j][i] = (float) w.get(j * patchSize + i, 0);
        }
    }

    imp1.setProcessor("test", new FloatProcessor(imageW));
    imp1.show();
    imp1.getCanvas().zoomIn(0, 0);
    imp1.getCanvas().zoomIn(0, 0);
    imp1.getCanvas().zoomIn(0, 0);
    imp1.getCanvas().zoomIn(0, 0);
    imp1.getCanvas().zoomIn(0, 0);
    imp1.getCanvas().zoomIn(0, 0);
    imp1.getCanvas().zoomIn(0, 0);

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