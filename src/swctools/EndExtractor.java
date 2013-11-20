package swctools;

import aux.ReadSWC;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 11/6/13
 * Time: 4:17 PM
 */
public class EndExtractor implements PlugIn {

	/*
        terminal call
        java -cp ~/critpoint/critpoint_.jar:ij.jar swctools.EndExtractor path-to-swc
     */
	public static void main(String args[]){

		/*
            check argument nr.
         */
		if (args.length!=1) {
			System.out.println("# GENERATE ENDPOINTS FROM SWC #");
			System.out.println("usage: path_to_SWC ");
			return;
		}

		/*
            check if the swc input file exists
         */
		String 	pathInSwc 	= new File(args[0]).getAbsolutePath();
		File 	fileSWC 	= new File(pathInSwc);

		if (!fileSWC.exists()) {
			System.out.println("file "+pathInSwc+" does not exist!");
			return;
		}

		/*
         set paths to outputs, same folder, keep the name with prefixes added
         */
		String parentDir = fileSWC.getParent() + File.separator;
		String name = fileSWC.getName().substring(0, fileSWC.getName().length()-4);
		// new names for outputs
		String pathOutEnd = parentDir + "END_" + name + ".swc";

		ReadSWC reader = new ReadSWC(pathInSwc);
		reader.exportEndpoints(pathOutEnd);

	}

	/*
        IJ call
     */
	public void run(String s) {

		GenericDialog gd = new GenericDialog("EXTRACT ENDPOINTS FROM RECONSTRUCTION");
		gd.addStringField("swc:", "path-to-swc", 50);
		gd.showDialog();
		if (gd.wasCanceled()) return;
		String pathInSwc = gd.getNextString();

		/*
            check if the swc input file exists
         */
		pathInSwc 	= new File(pathInSwc).getAbsolutePath();
		File 	fileSWC 	= new File(pathInSwc);

		if (!fileSWC.exists()) {
			System.out.println("file "+pathInSwc+" does not exist!");
			return;
		}

		/*
         set paths to outputs, same folder, keep the name with prefixes added
         */
		String parentDir = fileSWC.getParent() + File.separator;
		String name = fileSWC.getName().substring(0, fileSWC.getName().length()-4);
		// new names for outputs
		String pathOutEnd = parentDir + "END_" + name + ".swc";

		ReadSWC reader = new ReadSWC(pathInSwc);
		reader.exportEndpoints(pathOutEnd);

		IJ.log("exported \n" + pathOutEnd);

	}

}
