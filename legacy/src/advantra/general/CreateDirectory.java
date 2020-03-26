package advantra.general;

import java.io.*;

public class CreateDirectory {

	public static String createOneDir(String dir_name){
		
		File f = new File(dir_name);
		
		if(!f.exists()){
			
		try{
			// Create one directory
			
			boolean success = f.mkdir();
			if (!success) {
				System.out.println("Error: Directory: " + dir_name + "    .... was NOT created");
			}	
			
		}
			catch (Exception e){//Catch exception if any
				System.err.println("Error: " + e.getMessage());
		}
		
		}
		
		return f.getAbsolutePath();

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
