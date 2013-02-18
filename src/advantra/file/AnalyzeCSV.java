package advantra.file;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;

public class AnalyzeCSV {
	
	private String 	file_path;
	private int 	file_length;
	private int 	max_width = Integer.MIN_VALUE;
	private int 	min_width = Integer.MAX_VALUE;
	
	public AnalyzeCSV(String path){
		
		File f = new File(path);
		if(!f.exists()){
				System.err.println(""+f.getAbsolutePath()+" file does not exist!");
				return;
		}
		
		file_path		= f.getAbsolutePath();
		file_length		= 0;
		
		String read_line;
		
		// scan the file to check the length
		try {
			FileInputStream fstream 	= new FileInputStream(file_path);
			BufferedReader br 			= new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
			while ( (read_line = br.readLine()) != null ) {
				if(!read_line.trim().startsWith("#")) {
					String[] tokens = read_line.split( ",\\s*" );
					if(tokens.length>max_width) max_width = tokens.length;
					if(tokens.length<min_width) min_width = tokens.length;
					file_length++; // # are comments
				}
			}
			br.close();
			fstream.close();
		}
		catch (Exception e){
		  	System.err.println("Error: " + e.getMessage());
		}
	}
	
	/*
	 * reads first col_nr columns of every line
	 */
	public double[][] readLn(int col_nr){
		
		if(col_nr>min_width){
			return null;
		}
		
		double[][] cols = new double[file_length][col_nr];
		int read_line_number = 0;
		
		// scan the file
		try {
			FileInputStream fstream 	= new FileInputStream(file_path);
			BufferedReader br 			= new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
			String read_line;
			// should loop again the same amount of times it did at the line count
			while ((read_line = br.readLine()) != null) {
				
				if(!read_line.trim().startsWith("#")){
		    		
					String[] tokens = read_line.split( ",\\s*" );
		    		
					for (int i = 0; i < col_nr; i++) {
						cols[read_line_number][i] = Double.valueOf(tokens[i].trim()).doubleValue();
					}
					
					read_line_number++;
		    		
		    	}
			} // end looping the file
				
			br.close();
			fstream.close();	  		
			
		}
		catch (Exception e){
		  			System.err.println("Error: " + e.getMessage());
		}
		
		return cols;
	}

	/*
	 * read last column of every line
	 */
	public int[] readLastCol(){
		
		if(max_width!=min_width){
			System.out.println("Warning: Lines have different lengths! in file "+file_path);
		}
		
		int[] col = new int[file_length];
		int read_line_number = 0;
		
		// scan the file
		try {
			FileInputStream fstream 	= new FileInputStream(file_path);
			BufferedReader br 			= new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
			String read_line;
			// should loop again the same amount of times it did at the line count
			while ((read_line = br.readLine()) != null) {
				
				if(!read_line.trim().startsWith("#")){
		    		
					String[] tokens = read_line.split( ",\\s*" );
		    		
					col[read_line_number] = Integer.valueOf(tokens[tokens.length-1].trim()).intValue();
					
					read_line_number++;
		    		
		    	}
			} // end looping the file
				
			br.close();
			fstream.close();	  		
			
		}
		catch (Exception e){
		  			System.err.println("Error scanning the csv file: " + e.getMessage());
		}
		
		return col;
	}

	
	public int getLinesNr(){
		return file_length;
	}
	
	public int getMaxWidth(){
		return max_width;
	}
	
	public int getMinWidth(){
		return min_width;
	}

}