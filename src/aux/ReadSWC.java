package aux;

import java.io.*;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/23/13
 * Time: 5:02 PM
 */
public class ReadSWC {

    private String  swcFilePath = "";
    private int     fileLength  = 0;

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
            System.out.print("loading SWC...");
            while ( (read_line = br.readLine()) != null ) {
                if(!read_line.trim().startsWith("#")) {

                    fileLength++; // # are comments
                }
            }
            br.close();
            fstream.close();
        }
        catch (Exception e){
            System.err.println("Error: " + e.getMessage());
        }
        System.out.println(fileLength+" lines (nodes) found.");

    }

}
