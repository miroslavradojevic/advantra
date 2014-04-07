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
	public float minR = Float.POSITIVE_INFINITY, maxR = Float.NEGATIVE_INFINITY;
	public float minX = Float.POSITIVE_INFINITY, maxX = Float.NEGATIVE_INFINITY;
	public float minY = Float.POSITIVE_INFINITY, maxY = Float.NEGATIVE_INFINITY;
	public float minZ = Float.POSITIVE_INFINITY, maxZ = Float.NEGATIVE_INFINITY;

	// indexes
	public static int ID 		= 0;
	public static int TYPE 		= 1;
	public static int XCOORD 	= 2;
	public static int YCOORD 	= 3;
	public static int ZCOORD 	= 4;
	public static int RADIUS 	= 5;
	public static int MOTHER 	= 6;

    public ReadSWC(String swcFilePath) {

        swcFilePath = new File(swcFilePath).getAbsolutePath();

        if (!(new File(swcFilePath).exists())) {
            System.err.println(swcFilePath+" does not exist! class not initialized...");
            return;
        }

        this.swcFilePath = swcFilePath;
        this.fileLength = 0;

        try { // scan the file

			FileInputStream fstream 	= new FileInputStream(swcFilePath);
            BufferedReader br 			= new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
            String read_line;

			while ( (read_line = br.readLine()) != null ) {
                if(!read_line.trim().startsWith("#")) { // # are comments

                    fileLength++;

					// split values
					String[] 	readLn = 	read_line.trim().split("\\s+");

					float[] 	valsLn = 	new float[7]; // x, y, z, mother_index

                    valsLn[0] = Float.valueOf(readLn[ID].trim()).floatValue();  // id
                    valsLn[1] = Float.valueOf(readLn[TYPE].trim()).floatValue();  // type

					valsLn[2] = Float.valueOf(readLn[XCOORD].trim()).floatValue();  // x, y, z
					valsLn[3] = Float.valueOf(readLn[YCOORD].trim()).floatValue();
					valsLn[4] = Float.valueOf(readLn[ZCOORD].trim()).floatValue();

					valsLn[5] = Float.valueOf(readLn[RADIUS].trim()).floatValue();  // radius

					valsLn[6] = Float.valueOf(readLn[MOTHER].trim()).floatValue();  // mother idx

					nodes.add(valsLn);

					minR = (valsLn[RADIUS]<minR)? valsLn[RADIUS] : minR;
					maxR = (valsLn[RADIUS]>maxR)? valsLn[RADIUS] : maxR;

					minX = (valsLn[XCOORD]<minX)? valsLn[XCOORD] : minX;
					maxX = (valsLn[XCOORD]>maxX)? valsLn[XCOORD] : maxX;

					minY = (valsLn[YCOORD]<minY)? valsLn[YCOORD] : minY;
					maxY = (valsLn[YCOORD]>maxY)? valsLn[YCOORD] : maxY;

					minZ = (valsLn[ZCOORD]<minZ)? valsLn[ZCOORD] : minZ;
					maxZ = (valsLn[ZCOORD]>maxZ)? valsLn[ZCOORD] : maxZ;

                }
            }

			br.close();
            fstream.close();

		}
        catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
        }

    }

	public void print(){

		//		for (int ii=0; ii<nodes.size(); ii++) {
		//            for (int jj=0; jj<7; jj++) {
		//                System.out.print(nodes.get(ii)[jj]+"\t");
		//            }
		//            System.out.println();
		//		}

	}

	public void exportBifurcations(String swcBifExportPath){

		PrintWriter logWriter = null;

		try {
			logWriter = new PrintWriter(swcBifExportPath); logWriter.print(""); logWriter.close();
		} catch (FileNotFoundException ex) {}

		try {
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(swcBifExportPath, true)));
			//logWriter.println("# source "+ inSwc);
		} catch (IOException e) {}

		// will take loaded swc file and create a new one that extracts the bifurcation points
		int currId, currMotherId, laterMotherId, count;
		boolean isBif;

		// extraction
		for (int ii=0; ii<nodes.size(); ii++) {

			currId = Math.round(nodes.get(ii)[ID]);
			currMotherId = Math.round(nodes.get(ii)[MOTHER]);
			count = 0;
			isBif = false;

			// it is root, assign it as endpoint by default no need to loop
			if (currMotherId==-1) {
				for (int jj=ii+1; jj<nodes.size(); jj++) {

					laterMotherId = Math.round(nodes.get(jj)[MOTHER]);
					if (laterMotherId==currId) {
						count++;
						if (count==2) {
							isBif = true;
							break;
						}
					}

				}
			}

			if (isBif) {
				// add ii node to the output swc
				logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f -1",
													   currId,
													   5,  // type = 5
													   nodes.get(ii)[XCOORD],
													   nodes.get(ii)[YCOORD],
													   nodes.get(ii)[ZCOORD],
													   nodes.get(ii)[RADIUS]));

			}

		}

		logWriter.close();

		System.out.println(swcBifExportPath + " exported.");

	}

	public void exportEndpoints(String swcEndExportPath){

		PrintWriter logWriter = null;

		try {
			logWriter = new PrintWriter(swcEndExportPath);	logWriter.print(""); logWriter.close();
		} catch (FileNotFoundException ex) {}

		try {
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(swcEndExportPath, true)));
			//logWriter.println("# source "+ inSwc);
		} catch (IOException e) {}


		// will take current swc file and create a new one that marks the bifurcation points
		int currId, currMotherId, laterMotherId;
		boolean isEnd;

		// extraction
		for (int ii=0; ii<nodes.size(); ii++) {

			currId = Math.round(nodes.get(ii)[ID]);
			currMotherId = Math.round(nodes.get(ii)[MOTHER]);
			isEnd = true;

			if (currMotherId!=-1) {
				for (int jj=ii+1; jj<nodes.size(); jj++) {

					laterMotherId = Math.round(nodes.get(jj)[MOTHER]);
					if (laterMotherId==currId) {
						isEnd = false;
						break;
					}

				}
			}

			if (isEnd) {
				// add ii node to the output swc
				logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f -1",
													   currId,
													   6,  // type = 6
													   nodes.get(ii)[XCOORD],
													   nodes.get(ii)[YCOORD],
													   nodes.get(ii)[ZCOORD],
													   nodes.get(ii)[RADIUS]));

			}

		}

		logWriter.close();

		System.out.println(swcEndExportPath + " exported.");

	}

}
