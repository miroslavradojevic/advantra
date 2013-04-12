package advantra.plugins;

import java.io.File;

import javax.swing.JFileChooser;

import ij.IJ;
import ij.io.OpenDialog;

public class MyOpener {
	
	public static String open(String message) {
		
		IJ.showMessage(message);
		
		/*
		 * based on File_Opener.java ij source code 
		 * this sequence enables to browse any file 
		 * and get the path to it
		 */
		
		File dir=null;
		JFileChooser fc = null;
		try {fc = new JFileChooser();}
		catch (Throwable e) {return null;}
		if (dir==null) {
			String sdir = OpenDialog.getDefaultDirectory();
			if (sdir!=null)
				dir = new File(sdir);
		}
		if (dir!=null)
			fc.setCurrentDirectory(dir);
		int returnVal = fc.showOpenDialog(IJ.getInstance());
		if (returnVal!=JFileChooser.APPROVE_OPTION)
			return null;
		
		File train_dataset_file = fc.getSelectedFile();
		
		String train_dataset_path = fc.getSelectedFile().getAbsolutePath();
		
		if(!train_dataset_file.exists()){
			IJ.showMessage("file "+train_dataset_path+" does not exist");
			return null;
		}
		else{
			IJ.showMessage("opened    "+train_dataset_path);
		}
		
		return fc.getSelectedFile().getAbsolutePath();
		
	}

}
