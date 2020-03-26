package advantra.general;


import java.io.*;

public class DebugExport {
	
	private String  FILE_NAME;
	private FileWriter fw; 
	
	public DebugExport(String file_name){
		FILE_NAME = file_name;
		// open file
		try{fw = new FileWriter(FILE_NAME);} 
		catch(IOException exIO){System.out.printf("DebugExport:DebugExport():  \nCould not open DEBUG file for writing...");}
	}
	
	public void writeLine(String line_to_write){
		try
		{
			fw.write(line_to_write);
			fw.write("\n");
		} 
		catch(IOException exIO){System.out.printf("DebugExport:writeLine():  \nCould not write line to DEBUG file.");}
	}
	
	public void write(String contents_to_write){
		try
		{
			fw.write(contents_to_write);
		} 
		catch(IOException exIO){System.out.printf("DebugExport:write():  \nCould not write to DEBUG file.");}
	}
	
	public void writeMatrixMatlabForm(String variable_name, double[][] array_to_write){
		try
		{
			String array_string = variable_name + " = [ ... \n";
			for (int i = 0; i < array_to_write.length; i++) {
				for (int j = 0; j < array_to_write[0].length; j++) {
					array_string = array_string.concat(String.format("%.4f ", array_to_write[i][j]));
					if(j==array_to_write[0].length-1){	array_string = array_string.concat("; ");	}
				}
			}
			array_string = array_string.concat("];\n");
			fw.write(array_string);
		}
		catch(IOException exIO){System.out.printf("DebugExport:writeMatrixMatlabForm():  \nCould not write to DEBUG file.");}
	}
	
	/* VECTORS STORE */
	// DOUBLE
	/*
	public void writeDoubleVector(double[] input_vector_data, String input_vector_name, int input_vector_index){
		try{
			fw.write(input_vector_name+"{"+input_vector_index+",1} =... \n");
			fw.write("[");
			for(int i=0; i<input_vector_data.length; i++){
				fw.write(input_vector_data[i]+" ");
			}
			fw.write("];\n");
			
			} catch(IOException exIO){System.out.println("Couldn't write to debug file.");}
	}
	// INT
	public void writeIntVector(int[] input_vector_data, String input_vector_name, int input_vector_index){
		try{
			fw.write(input_vector_name+"{"+input_vector_index+",1} =... \n");
			fw.write("[");
			for(int i=0; i<input_vector_data.length; i++){
				fw.write(input_vector_data[i]+" ");
			}
			fw.write("];\n");
			
			} catch(IOException exIO){System.out.println("Couldn't write to debug file.");}
	}
	*/
	/* VALUES STORE */
	/*
	public void writeIntValue(int input_value, String value_name){
		try{
			fw.write(value_name+" =["+input_value+"];\n");
			} catch(IOException exIO){System.out.println("Couldn't write to debug file.");}
	}
	
	public void writeDoubleValue(double input_value, String value_name){
		try{
			fw.write(value_name+" =["+input_value+"];\n");
			} catch(IOException exIO){System.out.println("Couldn't write to debug file.");}
	}
	*/
	/* MATRICES STORE */
	/*
	// DOUBLE
	public void writeDoubleMatrix(double[][] input_matrix_data, String input_matrix_name, int input_matrix_index){
		try{
			fw.write(input_matrix_name+"{"+input_matrix_index+"} =... \n");
			fw.write("[");
			for(int i=0; i<input_matrix_data.length; i++){
				for(int j=0; j<input_matrix_data[0].length; j++){
					fw.write(input_matrix_data[i][j]+" ");
					if(j==input_matrix_data[0].length-1){fw.write("; ");}
				}
			}
			fw.write("];\n");
			} catch(IOException exIO){System.out.println("Couldn't write to debug file.");}
	}
	
	// INT
	public void writeIntMatrix(int[][] input_matrix_data, String input_matrix_name, int input_matrix_index){
		try{
			fw.write(input_matrix_name+"{"+input_matrix_index+"} =... \n");
			fw.write("[");
			for(int i=0; i<input_matrix_data.length; i++){
				for(int j=0; j<input_matrix_data[0].length; j++){
					fw.write(input_matrix_data[i][j]+" ");
					if(j==input_matrix_data[0].length-1){fw.write("; ");}
				}
			}
			fw.write("];\n");
			} catch(IOException exIO){System.out.println("Couldn't write to debug file.");}
	}
	*/
	
	public void closeDebug(){
			try{fw.close();} 
			catch(IOException exIO){System.out.printf("DebugExport:closeDebug():   \nCouldn't close debug file.");}
	}	
	
	
}
