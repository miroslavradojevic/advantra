import java.io.File;

/**
 * Created by miroslav on 18-4-15.
 */
public class SwcTest {

    public static void main(String[] args) {

        System.out.println("---");

//

        File f = new File(args[0]);
        if (!f.exists()) {System.out.println(args[0] + "\tis not existing file"); return;}

        if (!Toolbox.getFileExtension(f.getName()).equalsIgnoreCase("SWC")) {System.out.println(args[0] + "\tis not SWC"); return;}

        // read swc

        ReadSWC reader = new ReadSWC(f.getAbsolutePath());

//        System.out.println(f.getAbsolutePath());


    }

}
