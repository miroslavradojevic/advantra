 package advantra.general;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.ShortProcessor;

public class ImageConversions {
	
	/*
	 * different manipulations with ImageJ variables: 
	 * - 8-bit 
	 * - RGB
	 * - FloatImage
	 */

	// convert any ImagePlus to FloatImage
	public static ImagePlus ImagePlusToFloat(ImagePlus img){
		
		int im_width 	= img.getStack().getWidth();
		int im_height 	= img.getStack().getHeight();
		
		ImageStack im_stack = new ImageStack(im_width, im_height);
		
		for (int i = 1; i <= img.getStack().getSize(); i++) {
			// convert it first 
			// new ByteProcessor(img.getStack().getProcessor(i).createImage());
			// ImageProcessor fp = img.getStack().getProcessor(i).convertToFloat();
			// set FloatProcessor
			im_stack.addSlice(img.getStack().getProcessor(i).convertToFloat());
					
					
			//new FloatProcessor(im_width, im_height,img.getStack().getProcessor(i).getPixels()));
			
		}
		
		return new ImagePlus("Float_image", im_stack);
		
		
	}
	
	// convert any ImagePlus to rgb
	public static ImagePlus ImagePlusToRGB(ImagePlus img){
		
		ImageStack im_stack = new ImageStack(
				img.getStack().getWidth(), img.getStack().getHeight()); // container - layers iteratively added
		
		for (int i = 1; i <= img.getStack().getSize(); i++) {
			
			im_stack.addSlice(new ColorProcessor(img.getStack().getProcessor(i).createImage()));
			
		}
		
		return new ImagePlus("RGB_image", im_stack);
		
	}
	
	// any ImagePlus to Gray8
	
	// TODO: actually IJ.run(imageplus, "8-bit", ""); would do the same 
	public static ImagePlus ImagePlusToGray8(ImagePlus img){
		
		int w = img.getStack().getWidth();
		int h = img.getStack().getHeight();
		
		ImageStack im_stack = new ImageStack(w, h); // container - layers iteratively added
		
		byte[] array 	= new byte[w*h];
		byte[] red   	= new byte[w*h];
		byte[] green   	= new byte[w*h];
		byte[] blue   	= new byte[w*h];
		
		if(img.getType()==ImagePlus.COLOR_RGB){
			for (int i = 1; i <= img.getStack().getSize(); i++) {
				
				red 	= ((ColorProcessor)img.getStack().getProcessor(i)).getChannel(1);
				green 	= ((ColorProcessor)img.getStack().getProcessor(i)).getChannel(2);
				blue 	= ((ColorProcessor)img.getStack().getProcessor(i)).getChannel(3);
				
				for (int j = 0; j < array.length; j++) {
					array[j] = (byte)((int)Math.round(((int)(red[j]&0xff)+(int)(green[j]&0xff)+(int)(blue[j]&0xff))/3f)); // average
				}
				im_stack.addSlice(new ByteProcessor(img.getStack().getWidth(), img.getStack().getHeight(), array)); // img.getStack().getProcessor(i).createImage()
				
			}
			
			return new ImagePlus("Gray8_image", im_stack);
			
		}
		else{
			
			for (int i = 1; i <= img.getStack().getSize(); i++) {
				im_stack.addSlice(new ByteProcessor(img.getStack().getProcessor(i).createImage()));
			}
			
			return new ImagePlus("Gray8_image", im_stack);
			
		}
		
	}
	
	// convert it to int[][] format for regular processing, each 2d image is stored in a row-vector
	public static int[][] GraytoIntArray(ImagePlus img){
		
        int w = img.getStack().getWidth();
        int h = img.getStack().getHeight();
        int l = img.getStack().getSize();
        
        int[][]  img_stack = new int[l][w*h]; // where it will be stored
        
        // now check the type & read it
        if(img.getType()==ImagePlus.GRAY8){
        	
        	byte[] read_layer;
        	
        	for(int i = 1; i <= l; i++){
    			read_layer 	= (byte[])img.getStack().getPixels(i);
    			for (int j = 0; j < w*h; j++) {
    				img_stack[i-1][j] = read_layer[j] & 0xff;
    			}
    		}
        	
        }
        else if(img.getType()==ImagePlus.GRAY16){

        	short[] read_layer;
        	
        	for(int i = 1; i <= l; i++){
    			read_layer 	= (short[])img.getStack().getPixels(i);
    			for (int j = 0; j < w*h; j++) {
    				img_stack[i-1][j] = read_layer[j] & 0xffff;
    			}
    		}
        	
        }
        else{
        	
        	System.err.println("ImageConversions:toIntArray(): \nthis image type is not supported!");
			System.exit(1);	
        	
        }
		
		return img_stack;
	}

	public static int[][] GraytoIntArray(ImageStack img_stack_ij){
		
        int w = img_stack_ij.getWidth();
        int h = img_stack_ij.getHeight();
        int l = img_stack_ij.getSize();
        
        int[][]  img_stack_int = new int[l][w*h]; // where it will be stored
        
        // now check the type & read it
        if(img_stack_ij.getBitDepth()==8){
        	
        	byte[] read_layer;
        	
        	for(int i = 1; i <= l; i++){
    			read_layer 	= (byte[])img_stack_ij.getPixels(i);
    			for (int j = 0; j < w*h; j++) {
    				img_stack_int[i-1][j] = read_layer[j] & 0xff;
    			}
    		}
        	
        }
        else if(img_stack_ij.getBitDepth()==16){

        	short[] read_layer;
        	
        	for(int i = 1; i <= l; i++){
    			read_layer 	= (short[])img_stack_ij.getPixels(i);
    			for (int j = 0; j < w*h; j++) {
    				img_stack_int[i-1][j] = read_layer[j] & 0xffff;
    			}
    		}
        	
        }
        else{
        	
        	System.err.println("ImageConversions:toIntArray(): \nthis image type is not supported!");
			System.exit(1);	
        	
        }
		
		return img_stack_int;
	}
	
	// convert it to byte[][][] format for regular processing, each 2d image is stored in a row-vector
	// had to make it separate from gray because of rgb extraction (mainly getRGB() method that just couldn't combine r,g,b into one int)
	// byte[0][][] is red
	// byte[1][][] is green
	// byte[2][][] is blue
	public static byte[][][] RgbToByteArray(ImagePlus img){
		
        int w = img.getStack().getWidth();
        int h = img.getStack().getHeight();
        int l = img.getStack().getSize();
        
        byte[][][]  img_stack = new byte[3][l][w*h]; // where it will be stored
        
        // now check the type & read it
        if(img.getType()==ImagePlus.COLOR_RGB){
        	
        	for(int i = 1; i <= l; i++){
        		((ColorProcessor)(img.getStack().getProcessor(i))).getRGB(img_stack[0][i-1], img_stack[1][i-1], img_stack[2][i-1]);
    		}
        }
        else{
        	
        	System.err.println("ImageConversions:toIntArray(): \nother image type than rgb is not supported!");
			System.exit(1);	
        	
        }
		
		return img_stack;
	}

	// convert from int[][] to ByteImage (create new image)
	public static ImagePlus toGray8(int[][] img_stack, int width){
        
		int height 	= img_stack[0].length / width;
        int length 	= img_stack.length;
        
		ImageStack im_stack = new ImageStack(width, height); 	// container
		
		// this memory piece is feeding the ByteProcessor with data
		byte[][] values = new byte[length][width*height]; // will have to allocate the full size!!! just one 1d array is not enough
		
		// convert ints back to byte
		for (int i = 1; i <= length; i++) {
			
			for (int j = 0; j < width*height; j++) {
				
				values[i-1][j] = (byte)(img_stack[i-1][j]); // & 0xff
				
			}
			
			im_stack.addSlice(new ByteProcessor(width, height, values[i-1])); 
			
		}
		
		return new ImagePlus("Gray8Image", im_stack); // finally new image will be created with all the byte array memory assigned
		
	}
	
	public static ImagePlus toGray8(int[] img_stack, int width){
        
		int height 	= img_stack.length / width;
        
		ImageStack im_stack = new ImageStack(width, height); 	// container
		
		// this memory piece is feeding the ByteProcessor with data
		byte[] values = new byte[width*height]; // will have to allocate the full size!!! just one 1d array is not enough
		
		// convert ints back to byte

		for (int j = 0; j < width*height; j++) {
				
			values[j] = (byte)(img_stack[j]); 
				
		}
			
		im_stack.addSlice(new ByteProcessor(width, height, values)); 
		
		return new ImagePlus("Gray8Image", im_stack); // finally new image will be created with all the byte array memory assigned
		
	}	
	
	public static ImagePlus toGray8(byte[] img_stack, int width){
        
		int height 	= img_stack.length / width;
        
		ImageStack im_stack = new ImageStack(width, height); 	// container
		
		// this memory piece is feeding the ByteProcessor with data
		byte[] values = new byte[width*height]; // will have to allocate the full size!!! just one 1d array is not enough
		
		for (int j = 0; j < width*height; j++) {
				
			values[j] = img_stack[j]; 
				
		}
			
		im_stack.addSlice(new ByteProcessor(width, height, values)); 
		
		return new ImagePlus("Gray8Image", im_stack); // finally new image will be created with all the byte array memory assigned
		
	}

	public static ImagePlus toGray8(byte[][] values, int width){
        
		int height 	= values[0].length / width;
        int length 	= values.length;
        
		ImageStack im_stack = new ImageStack(width, height); 	// container
		
		// will have to allocate the full size!!! just one 1d array is not enough
		// this memory piece is feeding the ByteProcessor with data
		// byte[][] values = new byte[length][width*height];
		// convert it back to byte
		for (int i = 1; i <= length; i++) {
			
			im_stack.addSlice(new ByteProcessor(width, height, values[i-1])); 
		}
		
		return new ImagePlus("Gray8Image", im_stack); // finally new image will be created with all the byte array memory assigned
		
	}
	
	// convert from int[][] to ShortImage 
	public static ImagePlus toGray16(int[][] img_stack, int width){
        
		int h = img_stack[0].length / width;
        int l = img_stack.length;
        
		ImageStack im_stack = new ImageStack(width, h);  // container 
		
		short[][] read_layer = new short[l][width*h];
		
		// convert ints back to short
		for (int i = 1; i <= l; i++) {
			
			for (int j = 0; j < width*h; j++) {
				read_layer[i-1][j] = (short)(img_stack[i-1][j]); // & 0xff
			}
			
			im_stack.addSlice(new ShortProcessor(width, h, read_layer[i-1], java.awt.image.ColorModel.getRGBdefault()));
		}
		
		return new ImagePlus("Gray16Image", im_stack);
		
	}

	// convert from byte[][][] to ColorImage 
	public static ImagePlus toRgb(byte[][][] img_stack_array, int width){
        // this one is not really memory efficient
		int h = img_stack_array[0][0].length / width;
        //int l = img_stack_array[0].length;
        
		ImageStack im_stack = new ImageStack(width,h);  // container width, h, l
		System.out.println("!!!created ImageStack object");
		// java.awt.image.ColorModel.getRGBdefault() as argument?
		//int[] read_layer = new int[width*h];
		// convert int to int(rgb) - yes - it's a bit stupid! 
		// just to be consistent: img_stack and read_layer are the same thing right now
		
			for (int i = 1; i <= img_stack_array[0].length; i++) {
				im_stack.addSlice(String.format("slice%05d", i), new ColorProcessor(width, h), i-1);
				
				((ColorProcessor)im_stack.getProcessor(i)).setRGB(
						img_stack_array[0][i-1], 
						img_stack_array[1][i-1], 
						img_stack_array[2][i-1]
								);
				
			}
		
		ImagePlus out_img = new ImagePlus("AddedColorSpot", im_stack);
		System.out.println("!!!screated ImagePlus object");
		return out_img;
		
	}
	
	// sets (x,y,z) point of the image stack to a color expressed in jet colormap with given value (0-255)
	public static void setJet256ColorValue(ImagePlus img, int x, int y, int z, int value){
		
		x = (x<img.getHeight())?x:img.getHeight()-1;
		y = (y<img.getWidth() )?y:img.getWidth() -1;
		z = (z<img.getStack().getSize())?z:img.getStack().getSize()-1;

		x = (x>=0)?x:0;
		y = (y>=0)?y:0;
		z = (z>=0)?z:0;
		
		int w = img.getWidth();
		
		if((img.getType()!=ImagePlus.COLOR_RGB)){
			// set it to rgb
			img = ImageConversions.ImagePlusToRGB(img);
			System.out.println("converting ImagePlus to rgb...");
		}
		
		byte[] rgb_values = ColourTransf.Jet256(value);
		
		byte[][][] img_array = RgbToByteArray(img);
		
		img_array[0][z][ArrayHandling.sub2index_2d(x, y, w)] = rgb_values[0];
		img_array[1][z][ArrayHandling.sub2index_2d(x, y, w)] = rgb_values[1];
		img_array[2][z][ArrayHandling.sub2index_2d(x, y, w)] = rgb_values[2];
		
		for (int i = 1; i <= img.getStack().getSize(); i++) {
			
			((ColorProcessor)img.getStack().getProcessor(i)).setRGB(
					img_array[0][i-1], 
					img_array[1][i-1], 
					img_array[2][i-1]
							);
		}
	}
	
	// sets (x[],y[],z[]) points of the image stack to a color expressed in jet colormap with given value (0-255)
	public static void setJet256ColorValues(ImagePlus img, int[] x, int[] y, int[] z, int value){
		
		int[] x_index = new int[x.length];
		for (int i = 0; i < x.length; i++) {
			x_index[i] = (x[i]<img.getStack().getHeight())?x[i]:img.getStack().getHeight()-1;
			x_index[i] = (x_index[i]>=0)?x_index[i]:0;			
		}
		
		int[] y_index = new int[y.length];
		for (int i = 0; i < y.length; i++) {
			y_index[i] = (y[i]<img.getStack().getWidth())?y[i]:img.getStack().getWidth()-1;
			y_index[i] = (y_index[i]>=0)?y_index[i]:0;
		}
		
		int[] z_index = new int[z.length];
		for (int i = 0; i < z.length; i++) {
			z_index[i] = (z[i]<img.getStack().getSize())?z[i]:img.getStack().getSize()-1;
			z_index[i] = (z_index[i]>=0)?z_index[i]:0;
		}		
		
		int w = img.getWidth();
		
		if((img.getType()!=ImagePlus.COLOR_RGB)){
			// set it to rgb
			img = ImageConversions.ImagePlusToRGB(img);
			System.out.println("converting ImagePlus to rgb...");
		}
		
		byte[] rgb_values = ColourTransf.Jet256(value);
		
		byte[][][] img_array = RgbToByteArray(img);
		
		for(int i = 0; i < x.length; i++) {
			
			img_array[0][z_index[i]][ArrayHandling.sub2index_2d(x_index[i], y_index[i], w)] = rgb_values[0];
			img_array[1][z_index[i]][ArrayHandling.sub2index_2d(x_index[i], y_index[i], w)] = rgb_values[1];
			img_array[2][z_index[i]][ArrayHandling.sub2index_2d(x_index[i], y_index[i], w)] = rgb_values[2];
			
		}

		
		for (int i = 1; i <= img.getStack().getSize(); i++) {
			
			((ColorProcessor)img.getStack().getProcessor(i)).setRGB(
					img_array[0][i-1], 
					img_array[1][i-1], 
					img_array[2][i-1]
							);
		}
	}

	// sets (x[],y[],z[]) points of the image stack to a color expressed in jet colormap with given value (0-255)
	public static void setJet256ColorValues(ImagePlus img, double[] x, double[] y, double[] z, int value){
		
		int[] x_index = new int[x.length];
		for (int i = 0; i < x.length; i++) {
			x_index[i] = ((int)x[i]<img.getStack().getHeight())?(int)x[i]:img.getStack().getHeight()-1;
			x_index[i] = (x_index[i]>=0)?x_index[i]:0;			
		}
		
		int[] y_index = new int[y.length];
		for (int i = 0; i < y.length; i++) {
			y_index[i] = ((int)y[i]<img.getStack().getWidth())?(int)y[i]:img.getStack().getWidth()-1;
			y_index[i] = (y_index[i]>=0)?y_index[i]:0;
		}
		
		int[] z_index = new int[z.length];
		for (int i = 0; i < z.length; i++) {
			z_index[i] = ((int)z[i]<img.getStack().getSize())?(int)z[i]:img.getStack().getSize()-1;
			z_index[i] = (z_index[i]>=0)?z_index[i]:0;
		}		
		
		int w = img.getWidth();
		
		if((img.getType()!=ImagePlus.COLOR_RGB)){
			// set it to rgb
			img = ImageConversions.ImagePlusToRGB(img);
			System.out.println("converting ImagePlus to rgb...");
		}
		
		byte[] rgb_values = ColourTransf.Jet256(value);
		
		byte[][][] img_array = RgbToByteArray(img);
		
		for(int i = 0; i < x.length; i++) {
			
			img_array[0][z_index[i]][ArrayHandling.sub2index_2d(x_index[i], y_index[i], w)] = rgb_values[0];
			img_array[1][z_index[i]][ArrayHandling.sub2index_2d(x_index[i], y_index[i], w)] = rgb_values[1];
			img_array[2][z_index[i]][ArrayHandling.sub2index_2d(x_index[i], y_index[i], w)] = rgb_values[2];
			
		}

		
		for (int i = 1; i <= img.getStack().getSize(); i++) {
			
			((ColorProcessor)img.getStack().getProcessor(i)).setRGB(
					img_array[0][i-1], 
					img_array[1][i-1], 
					img_array[2][i-1]
							);
		}
	}

}
