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
	public ArrayList<float[]> nodes = new ArrayList<float[]>(); // 4x1 rows


    public ReadSWC(String swcFilePath) {

        swcFilePath = new File(swcFilePath).getAbsolutePath();

        if (!(new File(swcFilePath).exists())) {
            System.err.println(swcFilePath+" does not exist!");
            return;
        }

        this.swcFilePath = swcFilePath;
        this.fileLength = 0;

        // scan the file
        try {

			FileInputStream fstream 	= new FileInputStream(swcFilePath);
            BufferedReader br 			= new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
            String read_line;
            // check the length of the reconstruction first
//            System.out.print("loading SWC...");

			while ( (read_line = br.readLine()) != null ) {
                if(!read_line.trim().startsWith("#")) { // # are comments

                    fileLength++;

					// split values
					String[] 	readLn = 	read_line.trim().split("\\s+");


					//System.out.println("line:         "+read_line);
					//System.out.println("trimmed line: "+read_line.trim());
					//System.out.println("elements:     "+readLn.length);


					float[] 	valsLn = 	new float[5]; // x, y, z, mother_index

					valsLn[0] = Float.valueOf(readLn[2].trim()).floatValue();
					valsLn[1] = Float.valueOf(readLn[3].trim()).floatValue();
					valsLn[2] = Float.valueOf(readLn[4].trim()).floatValue();
					valsLn[3] = Float.valueOf(readLn[5].trim()).floatValue();

					valsLn[4] = Float.valueOf(readLn[6].trim()).floatValue();
					valsLn[4] = (valsLn[4]!=-1)?  valsLn[4]-1 : valsLn[4]  ; // decrease them so that they refer to index instead of the line

					nodes.add(valsLn);

                }
            }

			br.close();
            fstream.close();

		}
        catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
        }
//        System.out.println(nodes.size()+" lines (nodes) found.");

//		for (int ii=0; ii<nodes.size(); ii++) {
//			System.out.print(nodes.get(ii)[4]+", ");
//		}

    }

}
