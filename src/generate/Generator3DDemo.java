package generate;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 10/23/13
 * Time: 4:49 PM
 */
public class Generator3DDemo {

    public static void main(String args[]){

        float   k;
        float   snr;
        File    fileSWC;
        String  pathSWC;

        if (args.length!=3) {
            System.out.println("# GENERATE SYNTHETIC IMAGE FROM SWC #");
            System.out.println("usage: pathSWC k SNR");
            return;
        }

        pathSWC = new File(args[0]).getAbsolutePath();
        fileSWC = new File(pathSWC);

        if (!fileSWC.exists()) {
            System.out.println("file "+pathSWC+" does not exist!");
            return;
        }

        k       = Float.valueOf(args[1]);
        snr     = Float.valueOf(args[2]);

        System.out.println("loaded: \n"+pathSWC+"\n"+k+"\n"+snr);
        
        Generator3D g3D = new Generator3D();
        g3D.runFromSWC(pathSWC, k, snr);
        

    }

}
