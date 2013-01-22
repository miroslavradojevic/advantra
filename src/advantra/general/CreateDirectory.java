package advantra.general;

import java.io.*;

public class CreateDirectory {

	public static void createOneDir(String dir_name){
		try{
			// Create one directory
			boolean success = (new File(dir_name)).mkdir();
			if (success) {
				System.out.println("Directory: " + dir_name + "    ....created");
			}		
		}
			catch (Exception e){//Catch exception if any
				System.err.println("Error: " + e.getMessage());
		}
	}
	

	public static void createManyDirs(String dirs_names){
		try{
			// Create multiple directories
			boolean success = (new File(dirs_names)).mkdirs();
			if (success) {
				System.out.println("Directories: " + dirs_names + " created");
			}
		}
			catch (Exception e){//Catch exception if any
				System.err.println("Error: " + e.getMessage());
		}
	}

}
