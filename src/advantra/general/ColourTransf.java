package advantra.general;

public class ColourTransf {

	/*
	public static int getRed(int pix){
		return (int)(pix & 0xff0000)>>16;
	}
	
	public static int getGreen(int pix){
		return (int)(pix & 0x00ff00)>>8;
	}

	public static int getBlue(int pix){
		return (int)(pix & 0x0000ff);
	}
	
	public static int toPix(int red, int green, int blue){
		return ((red & 0xff) << 16)+((green & 0xff) << 8)+(blue & 0xff);
	}
	
	// conversion with arrays of values
	public static void getReds(int[] pix, int[] reds){
		
		if(pix.length != reds.length){
			System.err.println("Error in ColourTransf::getReds(): \n both procedure arguments have to have same length!");
			System.exit(1);
		}
		
		for (int i = 0; i < pix.length; i++) {
			reds[i] = (int)(pix[i] & 0xff0000)>>16;
		}
		
	}
	
	public static void getGreens(int[] pix, int[] greens){
		
		if(pix.length != greens.length){
			System.err.println("Error in ColourTransf::getGreens(): \n both procedure arguments have to have same length!");
			System.exit(1);
		}
		
		for (int i = 0; i < pix.length; i++) {
			greens[i] = (int)(pix[i] & 0x00ff00)>>8;
		}
		
	}

	public static void getBlues(int[] pix, int[] blues){
		
		if(pix.length != blues.length){
			System.err.println("Error in ColourTransf::getBlues(): \n both procedure arguments have to have same length!");
			System.exit(1);
		}
		
		for (int i = 0; i < pix.length; i++) {
			blues[i] = (int)(pix[i] & 0x0000ff);
		}
		
	}
	
	public static void toPixs(int[] reds, int[] greens, int[] blues, int[] pixs){
		if(
				(reds.length == greens.length) 	&& 
				(blues.length == pixs.length) 	&&
				(greens.length == blues.length)
				){
			for (int i = 0; i < pixs.length; i++) {
				pixs[i] = ((reds[i] & 0xff) << 16)+((greens[i] & 0xff) << 8)+(blues[i] & 0xff);
			}
		}
		else{
			System.err.println("Error in ColourTransf::toPixs(): \n all procedure arguments have to have the same length!");
			System.exit(1);
		}
			
		
	}
	*/
	
	public static byte[] Red(){
		
		byte[] rgb_values = new byte[3];
		
		rgb_values[0] = (byte)255;
	    rgb_values[1] = (byte)0;		 
	    rgb_values[2] = (byte)0;
	    
		return rgb_values;
		
	}
	
	public static byte[] Green(){
		
		byte[] rgb_values = new byte[3];
		
		rgb_values[0] = (byte)0;
	    rgb_values[1] = (byte)255;		 
	    rgb_values[2] = (byte)0;
	    
		return rgb_values;
		
	}
	
	public static byte[] Blue(){
		
		byte[] rgb_values = new byte[3];
		
		rgb_values[0] = (byte)0;
	    rgb_values[1] = (byte)0;		 
	    rgb_values[2] = (byte)255;
	    
		return rgb_values;
		
	}
	
	public static byte[] Jet256(int color_intensity){
		
		byte[] rgb_values = new byte[3];
		
		color_intensity = (color_intensity>255)?255:color_intensity;
		color_intensity = (color_intensity<0  )?0  :color_intensity;
		
		int[][] jet256 =
	    	 { 
	    		 { 0 , 0, 131 },
	    		 { 0 , 0, 135 },
	    		 { 0 , 0, 139 },
	    		 { 0 , 0, 143 },
	    		 { 0 , 0, 147 },
	    		 { 0 , 0, 151 },
	    		 { 0 , 0, 155 },
	    		 { 0 , 0, 159 },
	    		 { 0 , 0, 163 },
	    		 { 0 , 0, 167 },
	    		 { 0 , 0, 171 },
	    		 { 0 , 0, 175 },
	    		 { 0 , 0, 179 },
	    		 { 0 , 0, 183 },
	    		 { 0 , 0, 187 },
	    		 { 0 , 0, 191 },
	    		 { 0 , 0, 195 },
	    		 { 0 , 0, 199 },
	    		 { 0 , 0, 203 },
	    		 { 0 , 0, 207 },
	    		 { 0 , 0, 211 },
	    		 { 0 , 0, 215 },
	    		 { 0 , 0, 219 },
	    		 { 0 , 0, 223 },
	    		 { 0 , 0, 227 },
	    		 { 0 , 0, 231 },
	    		 { 0 , 0, 235 },
	    		 { 0 , 0, 239 },
	    		 { 0 , 0, 243 },
	    		 { 0 , 0, 247 },
	    		 { 0 , 0, 251 },
	    		 { 0 , 0, 255 },
	    		 { 0 , 4, 255 },
	    		 { 0 , 8, 255 },
	    		 { 0 , 12, 255 },
	    		 { 0 , 16, 255 },
	    		 { 0 , 20, 255 },
	    		 { 0 , 24, 255 },
	    		 { 0 , 28, 255 },
	    		 { 0 , 32, 255 },
	    		 { 0 , 36, 255 },
	    		 { 0 , 40, 255 },
	    		 { 0 , 44, 255 },
	    		 { 0 , 48, 255 },
	    		 { 0 , 52, 255 },
	    		 { 0 , 56, 255 },
	    		 { 0 , 60, 255 },
	    		 { 0 , 64, 255 },
	    		 { 0 , 68, 255 },
	    		 { 0 , 72, 255 },
	    		 { 0 , 76, 255 },
	    		 { 0 , 80, 255 },
	    		 { 0 , 84, 255 },
	    		 { 0 , 88, 255 },
	    		 { 0 , 92, 255 },
	    		 { 0 , 96, 255 },
	    		 { 0 , 100, 255 },
	    		 { 0 , 104, 255 },
	    		 { 0 , 108, 255 },
	    		 { 0 , 112, 255 },
	    		 { 0 , 116, 255 },
	    		 { 0 , 120, 255 },
	    		 { 0 , 124, 255 },
	    		 { 0 , 128, 255 },
	    		 { 0 , 131, 255 },
	    		 { 0 , 135, 255 },
	    		 { 0 , 139, 255 },
	    		 { 0 , 143, 255 },
	    		 { 0 , 147, 255 },
	    		 { 0 , 151, 255 },
	    		 { 0 , 155, 255 },
	    		 { 0 , 159, 255 },
	    		 { 0 , 163, 255 },
	    		 { 0 , 167, 255 },
	    		 { 0 , 171, 255 },
	    		 { 0 , 175, 255 },
	    		 { 0 , 179, 255 },
	    		 { 0 , 183, 255 },
	    		 { 0 , 187, 255 },
	    		 { 0 , 191, 255 },
	    		 { 0 , 195, 255 },
	    		 { 0 , 199, 255 },
	    		 { 0 , 203, 255 },
	    		 { 0 , 207, 255 },
	    		 { 0 , 211, 255 },
	    		 { 0 , 215, 255 },
	    		 { 0 , 219, 255 },
	    		 { 0 , 223, 255 },
	    		 { 0 , 227, 255 },
	    		 { 0 , 231, 255 },
	    		 { 0 , 235, 255 },
	    		 { 0 , 239, 255 },
	    		 { 0 , 243, 255 },
	    		 { 0 , 247, 255 },
	    		 { 0 , 251, 255 },
	    		 { 0 , 255, 255 },
	    		 { 4 , 255, 251 },
	    		 { 8 , 255, 247 },
	    		 { 12 , 255, 243 },
	    		 { 16 , 255, 239 },
	    		 { 20 , 255, 235 },
	    		 { 24 , 255, 231 },
	    		 { 28 , 255, 227 },
	    		 { 32 , 255, 223 },
	    		 { 36 , 255, 219 },
	    		 { 40 , 255, 215 },
	    		 { 44 , 255, 211 },
	    		 { 48 , 255, 207 },
	    		 { 52 , 255, 203 },
	    		 { 56 , 255, 199 },
	    		 { 60 , 255, 195 },
	    		 { 64 , 255, 191 },
	    		 { 68 , 255, 187 },
	    		 { 72 , 255, 183 },
	    		 { 76 , 255, 179 },
	    		 { 80 , 255, 175 },
	    		 { 84 , 255, 171 },
	    		 { 88 , 255, 167 },
	    		 { 92 , 255, 163 },
	    		 { 96 , 255, 159 },
	    		 { 100 , 255, 155 },
	    		 { 104 , 255, 151 },
	    		 { 108 , 255, 147 },
	    		 { 112 , 255, 143 },
	    		 { 116 , 255, 139 },
	    		 { 120 , 255, 135 },
	    		 { 124 , 255, 131 },
	    		 { 128 , 255, 128 },
	    		 { 131 , 255, 124 },
	    		 { 135 , 255, 120 },
	    		 { 139 , 255, 116 },
	    		 { 143 , 255, 112 },
	    		 { 147 , 255, 108 },
	    		 { 151 , 255, 104 },
	    		 { 155 , 255, 100 },
	    		 { 159 , 255, 96 },
	    		 { 163 , 255, 92 },
	    		 { 167 , 255, 88 },
	    		 { 171 , 255, 84 },
	    		 { 175 , 255, 80 },
	    		 { 179 , 255, 76 },
	    		 { 183 , 255, 72 },
	    		 { 187 , 255, 68 },
	    		 { 191 , 255, 64 },
	    		 { 195 , 255, 60 },
	    		 { 199 , 255, 56 },
	    		 { 203 , 255, 52 },
	    		 { 207 , 255, 48 },
	    		 { 211 , 255, 44 },
	    		 { 215 , 255, 40 },
	    		 { 219 , 255, 36 },
	    		 { 223 , 255, 32 },
	    		 { 227 , 255, 28 },
	    		 { 231 , 255, 24 },
	    		 { 235 , 255, 20 },
	    		 { 239 , 255, 16 },
	    		 { 243 , 255, 12 },
	    		 { 247 , 255, 8 },
	    		 { 251 , 255, 4 },
	    		 { 255 , 255, 0 },
	    		 { 255 , 251, 0 },
	    		 { 255 , 247, 0 },
	    		 { 255 , 243, 0 },
	    		 { 255 , 239, 0 },
	    		 { 255 , 235, 0 },
	    		 { 255 , 231, 0 },
	    		 { 255 , 227, 0 },
	    		 { 255 , 223, 0 },
	    		 { 255 , 219, 0 },
	    		 { 255 , 215, 0 },
	    		 { 255 , 211, 0 },
	    		 { 255 , 207, 0 },
	    		 { 255 , 203, 0 },
	    		 { 255 , 199, 0 },
	    		 { 255 , 195, 0 },
	    		 { 255 , 191, 0 },
	    		 { 255 , 187, 0 },
	    		 { 255 , 183, 0 },
	    		 { 255 , 179, 0 },
	    		 { 255 , 175, 0 },
	    		 { 255 , 171, 0 },
	    		 { 255 , 167, 0 },
	    		 { 255 , 163, 0 },
	    		 { 255 , 159, 0 },
	    		 { 255 , 155, 0 },
	    		 { 255 , 151, 0 },
	    		 { 255 , 147, 0 },
	    		 { 255 , 143, 0 },
	    		 { 255 , 139, 0 },
	    		 { 255 , 135, 0 },
	    		 { 255 , 131, 0 },
	    		 { 255 , 128, 0 },
	    		 { 255 , 124, 0 },
	    		 { 255 , 120, 0 },
	    		 { 255 , 116, 0 },
	    		 { 255 , 112, 0 },
	    		 { 255 , 108, 0 },
	    		 { 255 , 104, 0 },
	    		 { 255 , 100, 0 },
	    		 { 255 , 96, 0 },
	    		 { 255 , 92, 0 },
	    		 { 255 , 88, 0 },
	    		 { 255 , 84, 0 },
	    		 { 255 , 80, 0 },
	    		 { 255 , 76, 0 },
	    		 { 255 , 72, 0 },
	    		 { 255 , 68, 0 },
	    		 { 255 , 64, 0 },
	    		 { 255 , 60, 0 },
	    		 { 255 , 56, 0 },
	    		 { 255 , 52, 0 },
	    		 { 255 , 48, 0 },
	    		 { 255 , 44, 0 },
	    		 { 255 , 40, 0 },
	    		 { 255 , 36, 0 },
	    		 { 255 , 32, 0 },
	    		 { 255 , 28, 0 },
	    		 { 255 , 24, 0 },
	    		 { 255 , 20, 0 },
	    		 { 255 , 16, 0 },
	    		 { 255 , 12, 0 },
	    		 { 255 , 8, 0 },
	    		 { 255 , 4, 0 },
	    		 { 255 , 0, 0 },
	    		 { 251 , 0, 0 },
	    		 { 247 , 0, 0 },
	    		 { 243 , 0, 0 },
	    		 { 239 , 0, 0 },
	    		 { 235 , 0, 0 },
	    		 { 231 , 0, 0 },
	    		 { 227 , 0, 0 },
	    		 { 223 , 0, 0 },
	    		 { 219 , 0, 0 },
	    		 { 215 , 0, 0 },
	    		 { 211 , 0, 0 },
	    		 { 207 , 0, 0 },
	    		 { 203 , 0, 0 },
	    		 { 199 , 0, 0 },
	    		 { 195 , 0, 0 },
	    		 { 191 , 0, 0 },
	    		 { 187 , 0, 0 },
	    		 { 183 , 0, 0 },
	    		 { 179 , 0, 0 },
	    		 { 175 , 0, 0 },
	    		 { 171 , 0, 0 },
	    		 { 167 , 0, 0 },
	    		 { 163 , 0, 0 },
	    		 { 159 , 0, 0 },
	    		 { 155 , 0, 0 },
	    		 { 151 , 0, 0 },
	    		 { 147 , 0, 0 },
	    		 { 143 , 0, 0 },
	    		 { 139 , 0, 0 },
	    		 { 135 , 0, 0 },
	    		 { 131 , 0, 0 },
	    		 { 128 , 0, 0 }
	    		 };		
	    
	    rgb_values[0] = (byte)jet256[color_intensity][0];
	    rgb_values[1] = (byte)jet256[color_intensity][1];		 
	    rgb_values[2] = (byte)jet256[color_intensity][2];
	    
		return rgb_values;
		
	}
	
}
