package advantra.plugins;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;

import advantra.shapes.ConeCutoff;
import advantra.shapes.Sphere;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.NewImage;
import ij.io.FileSaver;
import ij.io.OpenDialog;
import ij.io.Opener;
import ij.plugin.PlugIn;
import ij.plugin.filter.GaussianBlur;
import ij.process.ImageProcessor;
import imagescience.random.PoissonGenerator;
import imagescience.utility.Progressor;

public class SynthetizeNeuron implements PlugIn {
	
	/*
	 * THIS PLUGIN WILL CREATE ARTIFICIAL NEURON IMAGE STACK
	 * USING swc FILE OF THE GIVEN RECONSTRUCTION
	 * 
	 * ORIGINAL PAPER
	 * xxxxxxxxxxxxxxxxxxxxxxxxxx
	 * 
	 * MATLAB IMPLEMENTATION AND FURTHER DETAILS AT:
	 * xxxxxxxxxxxxxxxxxxxxxxxxxx
	 * 
	 * INPUT IS PATH TO THE swc FILE
	 * FINAL OUTPUT IS THE BYTE NEURON STACK WITH 
	 * GAUSSIAN AND POISSON NOISE ADDED 
	 * 
	 * DEPENDING ON THE PLUGIN CALL IMAGE WITH DETECTIONS CAN BE 
	 * 1. DISPLAYED, (IMAGEJ PLUGIN CALL)
	 * 2. SAVED TO THE WORKING FOLDER (IF JAVA TERMINAL CALL)
	 * 	java -cp 
	 * /home/miroslav/fiji/jars/ij.jar:
	 * /home/miroslav/fiji/jars/imagescience.jar:
	 * /home/miroslav/fiji/plugins/Advantra_.jar 
	 * advantra.feature.PC_Extract_ 
	 * /home/miroslav/lena.png 4 6 3 2.1 0.55 2 0.5 10 -1
	 * 
	 * OR 
	 * 3. JUST EXPORTED AS ARRAY WITH INFORMATION FOR FURTHER USAGE (method get------ of BranchingPoints class)
	 */

	private static final String plugInTitle = "Synthetize neuron v0.1";
	private static final int MAX_R      = 15;
	
	public void run(String arg){
		
        // check if appropriate version of ImageJ is installed
        if (IJ.versionLessThan("1.24t"))
            return;

        // check Java
        if (System.getProperty("java.version").compareTo("1.6.0") < 0) {
            IJ.error("This plugin has been developed and tested with Java, version 1.7 and higher.\n" + "Please upgrade your JVM.");
            return;
        }

        IJ.run("Close All");
		
        
        
        Opener op = new Opener();
        op.open();
        
        String fileDir 	= OpenDialog.getLastDirectory();
        String fileName = OpenDialog.getLastName();
        String filePath = fileDir+fileName;
        
        // default param values
        int gaussianBlurSigma = 2;
        byte backgroundLevel = 20;
        double snr = 10.0;
        
		System.out.println(filePath);
		
				// dialog to enter input values
				GenericDialog gd = new GenericDialog(plugInTitle, IJ.getInstance());

				gd.addMessage("Gaussian blur :");
				gd.addNumericField( "sigma :", gaussianBlurSigma, 0, 5, "" );
				
				//gd.addNumericField("Gaussian blur sigma", 2,  0);	// default value, number of zeros
				
				gd.addMessage("Background level (Poisson) :");
				gd.addNumericField("intensity :", backgroundLevel, 0, 5, "");  
				
				gd.addMessage("SNR (Poisson) :");
				gd.addNumericField("ratio :", snr, 2, 5, "" );
				
				gd.showDialog();
				if (gd.wasCanceled()) return;
				
				
				gaussianBlurSigma = (int)gd.getNextNumber();
				backgroundLevel = (byte)gd.getNextNumber();
				snr = (double)gd.getNextNumber();

		// some counters used
		int fileLength	= 0;
		int currentLine = 0;
		 
		boolean doPoisson 	= true;
		boolean doGaussian 	= false;

  		// will assume that values stored in SWC are in pixels, not microns (influenced by DIADEM)
  		double coordinate_X_MIN = Double.POSITIVE_INFINITY; 
  		double coordinate_X_MAX = Double.NEGATIVE_INFINITY;
  		double coordinate_Y_MIN = Double.POSITIVE_INFINITY;
  		double coordinate_Y_MAX = Double.NEGATIVE_INFINITY;
  		double coordinate_Z_MIN = Double.POSITIVE_INFINITY;
  		double coordinate_Z_MAX = Double.NEGATIVE_INFINITY;

  		try{
  		
  			FileInputStream fstream 	= new FileInputStream(filePath);
  			BufferedReader br 			= new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
  			String strLine;
  			
  			// check the length of the reconstruction first
	  		System.out.println("Reading SWC...");
	  		while ( (strLine = br.readLine()) != null ) {
	  			if(!strLine.trim().startsWith("#")) fileLength++;
	  		}
		  	System.out.println("number of lines read :  "+fileLength);
  			
	  		// data to be read 
	  		double[][] coordinates 	= new double[fileLength][3];
	  		int[] label 			= new int[fileLength];
	  		int[] type 				= new int[fileLength];
	  		double[] radius 		= new double[fileLength]; 
	  		int[] parent 			= new int[fileLength];	
	  		System.out.println("allocated space for data to be read... ");
	  		
	  		// "reset" to beginning of file (discard old buffered reader)
	  		fstream.getChannel().position(0);
	  		br = new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
	  		
	  		System.out.println("reset reader to the beginning of the file...");
	  		// scan the lines of the file from the beginning again
	  		while ((strLine = br.readLine()) != null) {
	  			strLine = strLine.trim();
	  			if(!strLine.startsWith("#")){
		    		String[] tokens = strLine.split("\\s+");
		    		
		    		if(tokens.length == 7){
		    			
		    			label[currentLine] 			= Integer.valueOf(tokens[0].trim()).intValue();
		    			type[currentLine] 			= Integer.valueOf(tokens[1].trim()).intValue();
		    			
		    			// find X extrema
		    			if((coordinates[currentLine][0] = Double.valueOf(tokens[2].trim()).doubleValue())<coordinate_X_MIN) 
		    				{coordinate_X_MIN = coordinates[currentLine][0];}
		    			if(coordinates[currentLine][0]>coordinate_X_MAX)
		    				{coordinate_X_MAX = coordinates[currentLine][0];}
		    			
		    			// find Y extrema
		    			if((coordinates[currentLine][1] = Double.valueOf(tokens[3].trim()).doubleValue())<coordinate_Y_MIN) 
		    				{coordinate_Y_MIN = coordinates[currentLine][1];}
		    			if(coordinates[currentLine][1]>coordinate_Y_MAX)
		    				{coordinate_Y_MAX = coordinates[currentLine][1];}			    			
		    			
		    			// find Z extrema
		    			if((coordinates[currentLine][2] = Double.valueOf(tokens[4].trim()).doubleValue())<coordinate_Z_MIN) 
		    				{coordinate_Z_MIN = coordinates[currentLine][2];}
		    			if(coordinates[currentLine][2]>coordinate_Z_MAX)
		    				{coordinate_Z_MAX = coordinates[currentLine][2];}
		    			
		    			
		    			radius[currentLine]			= Double.valueOf(tokens[5].trim()).doubleValue();
		    			parent[currentLine] 		= Integer.valueOf(tokens[6].trim()).intValue();
		    			
		    			currentLine++; // counts the lines that were read
		    			//System.out.println("line "+currentLine+" read...");
		    		}
		    	}
	  		} // end looping the file
	  		
	  		br.close();
		    fstream.close();	  		
  		
  		// all the values are read now
  		System.out.println("Done with reading...");
  		
  		// set the lower borders
  		int startX=0; if(coordinate_X_MIN<0){startX=(int)Math.round(coordinate_X_MIN);}
  		int startY=0; if(coordinate_Y_MIN<0){startY=(int)Math.round(coordinate_Y_MIN);}
  		int startZ=0; if(coordinate_Z_MIN<0){startZ=(int)Math.round(coordinate_Z_MIN);}
  		startX = startX-2; startY = startY-2; startZ = startZ-2;
  		// set the higher borders
  		int endX=0; if(coordinate_X_MAX>0){endX=(int)Math.round(coordinate_X_MAX);}
  		int endY=0; if(coordinate_Y_MAX>0){endY=(int)Math.round(coordinate_Y_MAX);}
  		int endZ=0; if(coordinate_Z_MAX>0){endZ=(int)Math.round(coordinate_Z_MAX);}
  		
  		// extend the borders a bit
  		endX = endX+2; endY = endY+2; endZ = endZ+2;
  		
  		// total length in pixels
  		int x_length = endX-startX+1; int y_length = endY-startY+1; int z_length = endZ-startZ+1;
  		
  		System.out.println("Generating... ");
  		
  		// allocate the byte array that holds the whole stack
  		byte[][] voxels = new byte[z_length][x_length*y_length];
  		
  		// allocate byte image from scratch
  		ImagePlus current_image = NewImage.createByteImage("neuron_stack", x_length, y_length, z_length, NewImage.FILL_BLACK);
		ImageStack current_stack = current_image.getStack();
  		ImageProcessor current_ip;
  		
  		// store voxels in the stack
  		int stack_width = current_stack.getWidth();
  		
  		for (int layer = 1; layer <= z_length; layer++) {
  			
  			current_ip = current_stack.getProcessor(layer);
  			voxels[layer-1] = (byte[])current_ip.getPixels(); // voxels refers to the pixels in image stack now
  			
		}
  		
  		for (int i = 0; i < voxels.length; i++) {
			for (int j = 0; j < voxels[0].length; j++) {
		  		if(doPoisson){
		  			//int value = voxels[i][j] & 0xff;
		  			//value += backgroundLevel;
		  			voxels[i][j] = backgroundLevel; //(byte)value;
		  		}
		  		else
		  			voxels[i][j] = (byte)backgroundLevel;
			}
		}
  		
  		/* 
  		 * define the signal level to set
  		 * those pixels from the stack that belong to the spheres to I_o, the rest stays at the background level
  		 * mean is set wrt to defined snr and background level using the relation I_o/ sqrt(I_o + bg)=snr
  		 * I_o is obtained from the quadratic equation solution
		 */
  		
		double I_o = snr/2 * (snr + Math.sqrt(Math.pow(snr, 2) + 4 * backgroundLevel));
		
		System.out.println("Signal levels is: "+I_o);
		
  		/* 
  		 * mark the points that belong to the neuron - set them to logical 1 which is I_o actually 
		 * I_o is already defined and noise will be added according to these values
		 */

		Sphere 		sphere 		= new Sphere(0, 0, 0, 0); 
  		ConeCutoff 	coneCut 	= new ConeCutoff(0, 0, 0, 0, 0, 0, 0, 0);
  	
		int[]	sphere_voxel_vals	= new int	[Sphere.numberOfVoxInSphere(MAX_R)];
		int[][] sphere_voxel_coord	= new int[3][Sphere.numberOfVoxInSphere(MAX_R)];
  		
  		for (int i = 0; i < currentLine; i++) {
  			System.out.println("Processing line "+i+" / "+(currentLine-1));
  			
  			// check each SWC entry
  			
  			int pointX = (int)Math.round(coordinates[i][0]-startX);// in indexes
  			int pointY = (int)Math.round(coordinates[i][1]-startY);
  			int pointZ = (int)Math.round(coordinates[i][2]-startZ);
  			
  			
  			// extract indexes of the points around the sphere centered around the point
  			sphere.setCenter(pointX, pointY, pointZ); sphere.setR(radius[i]);
  			
  			int cnt = sphere.extractVox(current_image, sphere_voxel_coord, sphere_voxel_vals);//   extractVoxel   (current_stack);
  			
  			System.out.println("cords: "+sphere_voxel_coord.length+" x "+sphere_voxel_coord[0].length);
  			System.out.println("vals: "+sphere_voxel_coord.length);
  			System.out.println("DEBUG LINE ! "+cnt+"   radius = "+radius[i]);
  			
  			for (int j = 0; j < cnt; j++) {
  				
  				if(doPoisson){
  					System.out.println("first..."+sphere_voxel_coord[0][j]+" | "+sphere_voxel_coord[1][j]+" | "+sphere_voxel_coord[2][j]+"  |");
  					voxels[sphere_voxel_coord[2][j]][sphere_voxel_coord[1][j]*stack_width+sphere_voxel_coord[0][j]] = (byte)(I_o + backgroundLevel);
  				}
  				else{
  					voxels[sphere_voxel_coord[2][j]][sphere_voxel_coord[1][j]*stack_width+sphere_voxel_coord[0][j]] = (byte)(I_o + backgroundLevel);
  				}
			}
  			
  			System.out.println("DONE FOR THIS");
  			
  			// connect with parent node(s)
  			// find who the parent node is
  			if(parent[i]!=(-1)){
  				// find parent centroid and radius
  				for (int j = 0; j < currentLine; j++) {
					if(label[j]==parent[i]){
						
						// extract space between centroids - cutoff conic shape will be fit in between
						// i is current node, j is the parent that's been found
			  			coneCut.setCenter(coordinates[j][0]-startX, coordinates[j][1]-startY, coordinates[j][2]-startZ);
			  			coneCut.setCenter1(coordinates[i][0]-startX, coordinates[i][1]-startY, coordinates[i][2]-startZ);
			  			coneCut.setR(radius[j], radius[i]);
			  			int[][] cone_voxel_set = coneCut.extractVoxels(current_stack);
			  			
			  			// add those points that were extracted
			  			for (int k = 0; k < cone_voxel_set.length; k++) {
			  				if(doPoisson){
			  					voxels[cone_voxel_set[k][2]][cone_voxel_set[k][1]*stack_width+cone_voxel_set[k][0]] = (byte)(I_o + backgroundLevel);
			  				}
			  				else{
			  					voxels[cone_voxel_set[k][2]][cone_voxel_set[k][1]*stack_width+cone_voxel_set[k][0]] = (byte)(I_o + backgroundLevel);
			  				}
						}
			  			
					}
				}
  			}
  		} // end looping through the lines
  		
  		if(doPoisson){
  			// add poisson noise, define noise generator
	  		PoissonGenerator poissonGen = new PoissonGenerator();
	  		
	  		System.out.println("Adding Poisson noise...");
	  		
			Progressor pgs = new Progressor(); // progress bar
			pgs.display(true);
			pgs.steps(x_length*y_length*z_length); 
			pgs.start();
			
	  		for (int i = 0; i < voxels.length; i++) {
				for (int j = 0; j < voxels[0].length; j++) {
			  		int value = voxels[i][j] & 0xff;
					voxels[i][j] = (byte)poissonGen.next(value);
					pgs.step();
				}
			}
	  		pgs.stop();
  		}
  		
  		
  		if(doGaussian)
  		{
  		// add gaussian smooth
			for (int layer = 1; layer <= z_length; layer++) {
				
				current_ip = current_stack.getProcessor(layer);
				(new GaussianBlur()).blurGaussian(current_ip.convertToFloat(), gaussianBlurSigma, gaussianBlurSigma, 0.01);
//				if(!){
//					System.out.println("Something went wrong when bluring...");
//				}
			}
  		}	
  		
  		System.out.println("Showing the image...");
  		current_image.show();
  		
  		FileSaver synth_saver 	= new FileSaver(current_image);
  		boolean savedTiff;
  		savedTiff = synth_saver.saveAsTiffStack(fileDir+fileName+".tif");
  		if(!savedTiff){
  			System.out.println("Couldn't save TIFF.");
  		}
  		else{
  			System.out.println("TIF exported with the same name.");
  		}
    
  		}// try block for the file reader wraps the run method
  		catch (Exception e){
  			
  			System.err.println("Error: " + e.getMessage());
  			
  		}
  		
	}
	
	public static void main(String[] args){
		
		System.out.println("main()...");
		
	}
	
}
