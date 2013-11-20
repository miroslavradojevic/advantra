package aux;

import java.io.*;
import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/23/13
 * Time: 5:02 PM
 */
public class ReadSWC {

	// path to swc file
    private String  swcFilePath = "";
    // counter
	private int     fileLength  = 0;
	// list with nodes
	public ArrayList<float[]> nodes = new ArrayList<float[]>(); // 1x7 rows (swc format)
	public float minR = Float.MAX_VALUE, maxR = Float.MIN_VALUE;

    public ReadSWC(String swcFilePath) {

        swcFilePath = new File(swcFilePath).getAbsolutePath();

        if (!(new File(swcFilePath).exists())) {
            System.err.println(swcFilePath+" does not exist! class not initialized...");
            return;
        }

        this.swcFilePath = swcFilePath;
        this.fileLength = 0;

        // scan the file
        try {

			FileInputStream fstream 	= new FileInputStream(swcFilePath);
            BufferedReader br 			= new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
            String read_line;

			while ( (read_line = br.readLine()) != null ) {
                if(!read_line.trim().startsWith("#")) { // # are comments

                    fileLength++;

					// split values
					String[] 	readLn = 	read_line.trim().split("\\s+");

					float[] 	valsLn = 	new float[7]; // x, y, z, mother_index

                    valsLn[0] = Float.valueOf(readLn[0].trim()).floatValue();  // kept for consistency
                    valsLn[1] = Float.valueOf(readLn[1].trim()).floatValue();

					valsLn[2] = Float.valueOf(readLn[2].trim()).floatValue();  // x, y, z
					valsLn[3] = Float.valueOf(readLn[3].trim()).floatValue();
					valsLn[4] = Float.valueOf(readLn[4].trim()).floatValue();

					valsLn[5] = Float.valueOf(readLn[5].trim()).floatValue();  // radius

					valsLn[6] = Float.valueOf(readLn[6].trim()).floatValue();  // mother idx

//					valsLn[6] = Float.valueOf(readLn[6].trim()).floatValue();

                    // decrease them so that they refer to index instead of the line #
                    valsLn[6] = (valsLn[6]!=-1)?  valsLn[6]-1 : valsLn[6];
                    valsLn[0] = valsLn[0]-1;

                    // this way :
                    // 1 ..... -1
                    // 2 ..... 1
                    // 3 ..... 2

                    // becomes :
                    // 0 ..... -1
                    // 1 ..... 0
                    // 2 ..... 1


					nodes.add(valsLn);

					minR = (valsLn[5]<minR)? minR=valsLn[5] : minR ;
					maxR = (valsLn[5]>maxR)? maxR=valsLn[5] : maxR ;

                }
            }

			br.close();
            fstream.close();

		}
        catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
        }

//		// print the output
//		for (int ii=0; ii<nodes.size(); ii++) {
//            for (int jj=0; jj<7; jj++) {
//                System.out.print(nodes.get(ii)[jj]+"\t");
//            }
//            System.out.println();
//		}

    }

	public void exportBifurcations(String swcBifExportPath){

		/*
            initialize writer
         */
		PrintWriter logWriter = null;

		/*
			empty swc bif file
		 */
		try {
			logWriter = new PrintWriter(swcBifExportPath);
			logWriter.print("");
			logWriter.close();
		} catch (FileNotFoundException ex) {}

		/*
            initialize bif swc file
         */
		try {
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(swcBifExportPath, true)));
			//logWriter.println("# source "+ inSwc);
		} catch (IOException e) {}

		// will take current swc file and create a new one that marks the bifurcation points
		for (int ii=1; ii<nodes.size(); ii++) {
			// check whether it is a bifurcation (at least two other rows refer to it)

			int count = 0;
			int id_curr = Math.round(nodes.get(ii)[0]);
			for (int jj=ii+1; jj<nodes.size(); jj++) {
				int id_loop = Math.round(nodes.get(jj)[6]);
				if (id_loop==id_curr) {
					count++;
				}
			}

			if (count>=2) {
				// add the line for bifurcation point  (T=5 in swc format)
				logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f -1", Math.round(nodes.get(ii)[0]+1), 5, nodes.get(ii)[2], nodes.get(ii)[3], nodes.get(ii)[4], nodes.get(ii)[5]));
			}
		}

		logWriter.close();

		System.out.println(swcBifExportPath + " exported.");

	}

	public void exportEndpoints(String swcEndExportPath){

		/*
            initialize writer
         */
		PrintWriter logWriter = null;

		/*
			empty swc bif file
		 */
		try {
			logWriter = new PrintWriter(swcEndExportPath);
			logWriter.print("");
			logWriter.close();
		} catch (FileNotFoundException ex) {}

		/*
            initialize bif swc file
         */
		try {
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(swcEndExportPath, true)));
			//logWriter.println("# source "+ inSwc);
		} catch (IOException e) {}

		// will take current swc file and create a new one that marks the bifurcation points
		for (int ii=1; ii<nodes.size(); ii++) {
			// check whether it is a endpoint (zero other rows refer to it)

			int count = 0;
			int id_curr = Math.round(nodes.get(ii)[0]);
			for (int jj=ii+1; jj<nodes.size(); jj++) {
				int id_loop = Math.round(nodes.get(jj)[6]);
				if (id_loop==id_curr) {
					count++;
				}
			}

			if (count==0) {
				// add the line for endpoint  (T=6 in swc format)
				logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f -1", Math.round(nodes.get(ii)[0]+1), 6, nodes.get(ii)[2], nodes.get(ii)[3], nodes.get(ii)[4], nodes.get(ii)[5]));
			}
		}

		logWriter.close();

		System.out.println(swcEndExportPath + " exported.");

	}

}
