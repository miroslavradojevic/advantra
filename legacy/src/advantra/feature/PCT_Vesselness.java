 package advantra.feature;

//import Jama.EigenvalueDecomposition;
// import Jama.Matrix;
import flanagan.math.Matrix;
import ij.ImagePlus;
import ij.gui.NewImage;
import ij.io.FileSaver;
import ij.process.ImageProcessor;
import imagescience.feature.Smoother;
import imagescience.image.FloatImage;


public class PCT_Vesselness {
	
	private ImageProcessor 	imp;
	private String 			im_name;
	private float[] 		PCT_Vess;
	private float[][] 		PCT_VectorField;
	// private boolean 		isDone;
	
	public PCT_Vesselness(ImageProcessor imp, String im_name){ // constructor
		
		this.imp 		= imp;
		this.im_name 	= im_name;
		
		this.PCT_Vess		= new float[imp.getWidth()*imp.getHeight()];
		this.PCT_VectorField= new float [2][imp.getWidth()*imp.getHeight()];
		
		//this.isDone		= false;
		
	}
	
	public void setImage(ImageProcessor imp, String im_name){// to change the input image if necessary
		
		this.imp 		= imp;
		this.im_name 	= im_name;
		
		this.PCT_Vess		= new float[imp.getWidth()*imp.getHeight()];
		this.PCT_VectorField= new float [2][imp.getWidth()*imp.getHeight()];
		
		// this.isDone		= false;
		
	}
	
	public float[] getVess(){
		
		return PCT_Vess;
		
	}
	
	public float[][] getVectField(){
		
		return PCT_VectorField;
		
	}
	
	
	public void run(String scales, int norients, double k, double b, double c, boolean isCalledFromImageJ, boolean saveIt){
		
		// extract the scales from the string argument scales
		String[] after_splitting = scales.split(",\\s*");
		
		// split using regex comma and 0 or any other number of spaces
		
		int number_of_scales = after_splitting.length;
		double[] sc = new double[number_of_scales];
		
		int im_height = imp.getHeight();
		int im_width  = imp.getWidth();
		
		if(number_of_scales>0){
			
			for (int i = 0; i < number_of_scales; i++) {
				
				sc[i] = Double.valueOf(after_splitting[i]);
			
			}
			
		}
		else{
			
			System.err.println("There are not enough scales.");
			// isDone = false; 
			return;
		
		}
		
		FloatImage im_scaled 	= new FloatImage(new ImagePlus("ScaledImage", 	imp.convertToFloat()));
		FloatImage im_orig 		= new FloatImage(new ImagePlus("OrigImage", 	imp.convertToFloat()));
		
		Smoother smooth = new Smoother();
		
		float[][][] PC = new float[number_of_scales][norients][im_height*im_width];
		
		// set params
		PC_Params pc_parameters = new PC_Params();
		pc_parameters.norient = norients;
		pc_parameters.k = (float) k;

		for (int i = 0; i < number_of_scales; i++) { 
			
			im_scaled = im_orig;
			smooth.gauss(im_scaled, sc[i]);
			
			String name = "_"+im_name+"_scale_"+Double.toString(sc[i]);
			
			// run(ImageProcessor imp, String im_title, PC_Params params, boolean isCalledFromImageJ, boolean saveIt)
			PC[i] = PC_Extract.run(im_scaled.imageplus().getProcessor(), name, pc_parameters, false, false);

		}
		
		// orientations n_o = [ -cos(theta_o) -sin(theta_o) ] 
		double[][] n_o = new double [2][norients];
		for (int ort = 0; ort < norients; ort++) {
			
			n_o[0][ort] = - Math.cos(ort*Math.PI/norients);//.setElement(0, orientation, );
			n_o[1][ort] = - Math.sin(ort*Math.PI/norients);//.setElement(1, orientation, );
			// this way the range is [-pi, 0) and eig vectors should be wrapped within that range
		}
		
		Matrix tensor 			= new Matrix(2, 2); // using flanagan's library
		Matrix tensor_incr		= new Matrix(2, 2);
		
		double [] 		eigVals = new double[2];
		double [][] 	eigVect = new double[2][2];
		double       	tempVal; // for swapping 
		
		float pct_vesselness; // to store intermediate result
		
		// EigenvalueDecomposition eig_dec = new EigenvalueDecomposition(tensor);
		// Matrix eig_dec = new Matrix(2,2);
		// eigen decomp. for each pixel location and each scale
		for (int sca = 0; sca < number_of_scales; sca++) {
			// extract vesselness for particular point
			for (int idx = 0; idx < im_height*im_width; idx++) {
				
				// vesselness
				tensor = tensor.times(0); // set it to zero before summing
					
				for (int ort = 0; ort < norients; ort++) {
					
					// n_o*n_o^T - unit(2,2) tensor calculation
					tensor_incr.setElement(0, 0, PC[sca][ort][idx] * (n_o[0][ort]*n_o[0][ort]-1));
					tensor_incr.setElement(0, 1, PC[sca][ort][idx] *  n_o[0][ort]*n_o[1][ort]);
					tensor_incr.setElement(1, 0, PC[sca][ort][idx] *  n_o[0][ort]*n_o[1][ort]);
					tensor_incr.setElement(1, 1, PC[sca][ort][idx] * (n_o[1][ort]*n_o[1][ort]-1));
						
					tensor.plusEquals(tensor_incr);
						
				}

				// eigen analysis
				//EigenvalueDecomposition eig_dec = new EigenvalueDecomposition(tensor);
				eigVals = tensor.getEigenValues();
				eigVect = tensor.getEigenVector();
				// eigVals 	= eig_dec.getRealEigenvalues();
				// eigVect		= eig_dec.getV().getArray(); // The columns of eigVect represent the eigenvectors
				
				if( Math.abs(eigVals[1]) < Math.abs(eigVals[0]) ){
					
					// swap values
					tempVal = eigVals[0];
					eigVals[0] = eigVals[1];
					eigVals[1] = tempVal;
					// swap vec columns
					tempVal = eigVect[0][0];
					eigVect[0][0] = eigVect[0][1];
					eigVect[1][0] = tempVal;
					//
					tempVal = eigVect[1][0];
					eigVect[1][0] = eigVect[1][1];
					eigVect[1][1] = tempVal;
					
				}
				if(eigVals[1]>0){
					
					pct_vesselness = 0;
				}
				else{
					
					pct_vesselness = 
								(float) ((float) Math.exp(-(eigVals[0]*eigVals[0])/(2*b*b*eigVals[1]*eigVals[1]))*
										(1-Math.exp(-(eigVals[0]*eigVals[0]+eigVals[1]*eigVals[1])/(2*c*c))));
				}
				
				if(pct_vesselness>PCT_Vess[idx]){
					
					PCT_Vess[idx] = pct_vesselness;
					PCT_VectorField[0][idx] = (float) eigVect[0][0];
					PCT_VectorField[1][idx] = (float) eigVect[1][0];
					
				}
			} // for idx
		} // for sc
		
		//isDone = true;
		
		//} // if(!isDone)
		
		if(isCalledFromImageJ || saveIt){
			
			ImagePlus out_image = NewImage.createFloatImage (
					"pctvess", imp.getWidth(), imp.getHeight(), 1, NewImage.FILL_BLACK);
			
			out_image.getProcessor().setPixels(PCT_Vess);
			
			if(isCalledFromImageJ){
				
				// show it
				out_image.show();
			
			}
			
			if(saveIt){
				
				// save it to HD... 
				String output_title		= 
						im_name+"_pct_vess_orts"+Integer.toString(norients)+"_k"+Double.toString(k)+"_b"+Double.toString(b)+"_c"+Double.toString(c)+".png";
				boolean savedTiff = (new FileSaver(out_image)).saveAsPng(output_title); 
				
				if(!savedTiff){
					
					System.out.println("Couldn't save the png.");
			  		
				}
			  	else{
			  			
			  		System.out.println(output_title+" exported.");
			  		
			  	}
			}
		}
	
	} // run()

}
