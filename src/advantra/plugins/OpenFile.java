package advantra.plugins;

import java.io.File;

import javax.swing.JFileChooser;

import ij.IJ;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

public class OpenFile  implements PlugIn {
	
	String 	train_dataset_path    = "";
	File 	train_dataset_file;
	
	public void run(String arg0) {
		
		IJ.showMessage("Open some file");
		
		/*
		 * based on File_Opener.java ij source code 
		 * this sequence enables to browse any file 
		 * and get the path to it
		 */
		
		File dir=null;
		JFileChooser fc = null;
		try {fc = new JFileChooser();}
		catch (Throwable e) {return;}
		if (dir==null) {
			String sdir = OpenDialog.getDefaultDirectory();
			if (sdir!=null)
				dir = new File(sdir);
		}
		if (dir!=null)
			fc.setCurrentDirectory(dir);
		int returnVal = fc.showOpenDialog(IJ.getInstance());
		if (returnVal!=JFileChooser.APPROVE_OPTION)
			return;
		train_dataset_file = fc.getSelectedFile();
		train_dataset_path = fc.getSelectedFile().getAbsolutePath();
		if(!train_dataset_file.exists()){
			IJ.showMessage("file "+train_dataset_path+" does not exist");
			return;
		}
		else{
			IJ.showMessage("opened file path"+train_dataset_path);
		}
		
	}

}
