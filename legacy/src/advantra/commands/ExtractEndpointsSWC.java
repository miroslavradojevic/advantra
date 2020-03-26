package advantra.commands;

import advantra.file.AnalyzeSWC;

public class ExtractEndpointsSWC {

	public static void main(String[] args){

		String swc_path = "";
		String img_path = "";
		
		switch(args.length){
		case 2:
			swc_path = args[0];
			img_path = args[1];
			break;
		default:
			System.err.println(
					"Extract bifurcations... they will be saved in a file with positions & neurite radius at the position \n" +
					"Additional image is exported visualizing what was saved. \n" +
					"1 - swc file path\n" +
					"2 - image stack path\n");
			System.exit(1);
		}
		
		AnalyzeSWC a = new AnalyzeSWC(swc_path);
		a.load();
		int visualization_color = 255;
		a.extractEndpoints(img_path, visualization_color);
 		
		System.out.println("done.");

	}
	
}
