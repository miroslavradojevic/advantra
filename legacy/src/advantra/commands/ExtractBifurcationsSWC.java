package advantra.commands;

import java.io.File;

import advantra.file.AnalyzeSWC;

public class ExtractBifurcationsSWC {

	public static void main(String[] args){

		String swc_path = "";
		String img_path = "";
		
		switch(args.length){
		case 2:
			swc_path = args[0];
			img_path = args[1];
			break;
		default:
			System.err.println("Extract bifurcations... they will be saved in a file with positions \n" +
					"1 - swc file path\n" +
					"2 - image stack path\n");
			System.exit(1);
		}
		
		// check
		swc_path = (new File(args[0])).getAbsolutePath();
		if(!(new File(swc_path)).exists()){
			System.err.println(swc_path+"  does not exist!");
			return;
		}
		
		img_path = (new File(args[1])).getAbsolutePath();
		if(!(new File(img_path)).exists()){
			System.err.println(img_path+"  does not exist!");
			return;
		}
		System.out.println("loaded files: \nSWC: "+swc_path+" ");
		System.out.println("TIF: "+img_path);
		
		AnalyzeSWC a = new AnalyzeSWC(swc_path);
		a.load();
		int visualization_color = 255;
		a.extractBifurcations(img_path, visualization_color); // extract both text file & visualize in an image stack
 		
		System.out.println("done.");

	}
	
}
