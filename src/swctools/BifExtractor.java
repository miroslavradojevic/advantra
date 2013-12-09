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
 * Time: 2:54 PM
 */
public class BifExtractor implements PlugIn {

	/*
        terminal call
        java -cp ~/critpoint/critpoint_.jar:ij.jar swctools.BifExtractor path-to-swc
     */
	public static void main(String args[]){

		/*
            check argument nr.
         */
		if (args.length!=1) {
			System.out.println("# GENERATE BIFURCATIONS FROM SWC #");
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
		String pathOutBif = parentDir + "BIF_" + name + ".swc";

		ReadSWC reader = new ReadSWC(pathInSwc);
		reader.exportBifurcations(pathOutBif);

	}

	/*
        IJ call
     */
	public void run(String s) {

		GenericDialog gd = new GenericDialog("EXTRACT BIFURCATIONS FROM RECONSTRUCTION");
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
		String pathOutBif = parentDir + "BIF_" + name + ".swc";

		ReadSWC reader = new ReadSWC(pathInSwc);
		reader.exportBifurcations(pathOutBif);

		IJ.log("exported \n" + pathOutBif);
		// TODO add the option to export them in csv list

	}

}
