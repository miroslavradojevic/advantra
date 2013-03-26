package advantra.critpoint;

import java.awt.Button;
import java.awt.FlowLayout;
import java.awt.Panel;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileFilter;

import MLdetection.FeaturePool;

import ij.IJ;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

public class TrainPatches implements PlugIn, ActionListener {

//	Button posButton;
	int nr_runs;
	
	public void run(String arg0) {
		
		nr_runs = (int)Prefs.get("advantra.critpoint.TrainPatches.nr_runs", 3);
		
		GenericDialog gd = new GenericDialog("Train Patches");
	    //dlg.addNumericField("Size of the scanning window : ", patchSize, 0);
	    gd.addNumericField("Number of adaboost runs : ", nr_runs, 0);
//	    dlg.addNumericField("Search Step for weak classifier : ", weakStep, 0);
	    //dlg.addPanel(makeButtonPanel_image(dlg));
	    gd.addPanel(makeButtonPanel_train());
//	    dlg.addPanel(makeButtonPanel_training_n(dlg));
	    gd.showDialog();
	    if (gd.wasCanceled()) {
	        return;
	    }

	   // patchSize = (int) dlg.getNextNumber();
	    nr_runs = (int) gd.getNextNumber();
//	    weakStep = (int) dlg.getNextNumber();
		
		// create the feature pool 
	    FeaturePool featurePool = new FeaturePool(5);
	    int fSize = featurePool.getSize();
	    IJ.log("the fearute pool of size " + fSize + " is created");
		
	}
	
	private Panel makeButtonPanel_train() { //GenericDialog gd
	    
		Panel buttons = new Panel();
	    buttons.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 1));

	    Button posButton;
	    posButton = new Button("Open stack with positive trainset");
	    posButton.addActionListener(this);
	    buttons.add(posButton);

//	    fnTrainP = new TextField(fileNamePositiveExamples, 80);
//	    fnTrainP.addActionListener(this);
//	    buttons.add(fnTrainP);

	    return buttons;
	}

	public void actionPerformed(ActionEvent e) {
		if(e.getSource()==posButton){
			JFileChooser fc = new JFileChooser();
			int returnVal = fc.showOpenDialog(null);
	        if (returnVal != JFileChooser.APPROVE_OPTION) {
	            return;
	        }
	        String fileName = fc.getSelectedFile().getAbsolutePath();
	        IJ.showStatus("Opened " + fileName);
		}
	}

}
