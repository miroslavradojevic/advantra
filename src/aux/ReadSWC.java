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
                    valsLn[2] = Float.valueOf(readLn[2].trim()).floatValue();
					valsLn[3] = Float.valueOf(readLn[3].trim()).floatValue();
					valsLn[4] = Float.valueOf(readLn[4].trim()).floatValue();
					valsLn[5] = Float.valueOf(readLn[5].trim()).floatValue();
					valsLn[6] = Float.valueOf(readLn[6].trim()).floatValue();

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

                }
            }

			br.close();
            fstream.close();

		}
        catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
        }

		for (int ii=0; ii<nodes.size(); ii++) {
            for (int jj=0; jj<7; jj++) {
                System.out.print(nodes.get(ii)[jj]+"\t");
            }
            System.out.println();
		}

    }

}
