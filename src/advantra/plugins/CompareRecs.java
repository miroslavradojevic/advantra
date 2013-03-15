package advantra.plugins;

import java.io.File;

import javax.swing.JFileChooser;

import advantra.file.AnalyzeCSV;
import ij.IJ;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

public class CompareRecs implements PlugIn {

	String 	path1, path2;// = "";
	File	file1, file2;// = "";
	
	public void run(String arg0) {
		
		// open file 1
		
		File dir=null;
		JFileChooser fc = null;
		try {fc = new JFileChooser();}
		catch (Throwable e) {IJ.error("This plugin requires Java 2 or Swing."); return;}
		//fc.setMultiSelectionEnabled(true);
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
		file1 = fc.getSelectedFile();
		path1 = fc.getSelectedFile().getAbsolutePath();
		if(!file1.exists()){
			IJ.showMessage("file "+path1+" does not exist");
			return;
		}
		else{
			IJ.showMessage("Chosen file path"+path1);
		}
		
		// open file 2
		
		dir=null;
		fc = null;
		try {fc = new JFileChooser();}
		catch (Throwable e) {IJ.error("This plugin requires Java 2 or Swing."); return;}
		//fc.setMultiSelectionEnabled(true);
		if (dir==null) {
			String sdir = OpenDialog.getDefaultDirectory();
			if (sdir!=null)
				dir = new File(sdir);
		}
		if (dir!=null)
			fc.setCurrentDirectory(dir);
		returnVal = fc.showOpenDialog(IJ.getInstance());
		if (returnVal!=JFileChooser.APPROVE_OPTION)
			return;
		file2 = fc.getSelectedFile();
		path2 = fc.getSelectedFile().getAbsolutePath();
		if(!file1.exists()){
			IJ.showMessage("file "+path2+" does not exist");
			return;
		}
		else{
			IJ.showMessage("Chosen file path"+path2);
		}
		
		IJ.showMessage("Compare files:\n"+path2+"\n and \n"+path1);
		
		AnalyzeCSV a1 = new AnalyzeCSV(path1);
		AnalyzeCSV a2 = new AnalyzeCSV(path2);
		
		System.out.println("done!");
				
		
	}
	
	

}
