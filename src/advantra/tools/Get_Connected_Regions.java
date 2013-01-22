package advantra.tools;

import java.awt.image.ColorModel;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Random;

import amira.AmiraParameters;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.NewImage;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.ImageCalculator;
import ij.process.ByteProcessor;

public class Get_Connected_Regions {
	
	boolean diagonal;
	boolean display;
	boolean showResults;
	boolean mustHaveSameValue;
	boolean startFromPointROI;
	boolean autoSubtract;
	double valuesOverDouble;
	double minimumPointsInRegionDouble;
	int stopAfterNumberOfRegions;
	
	boolean pleaseStop;
	
	ImagePlus img;
	
	public Get_Connected_Regions(ImagePlus img){
		
//		boolean diagonal 	= true;
//		boolean display 	= false;
//		boolean showResults = false;
//		boolean mustHaveSameValue = true;
//		boolean startFromPointROI = false;
//		boolean autoSubtract = false;
//		double 	valuesOverDouble = 0.0;
//		double 	minimumPointsInRegionDouble = 0.0;
//		int 	stopAfterNumberOfRegions = -1;
		
		this.img = img;
		
		if(this.img.getType()!=ImagePlus.GRAY8){
			IJ.error("The image must be either 8 bit.");
			return;
		}
		
//		boolean pleaseStop = false;
		
	}
	
	public void run(){
		
		int type = img.getType();
		
		boolean byteImage = true;
		//boolean startAtMaxValue = false; // !mustHaveSameValue;
		//ImageCalculator iCalc = new ImageCalculator();
		
		int point_roi_x = -1;
		int point_roi_y = -1;
		int point_roi_z = -1;

		int width 	= img.getWidth();
		int height 	= img.getHeight();
		int depth 	= img.getStackSize();

		if (width * height * depth > Integer.MAX_VALUE) {
			IJ.error("This stack is too large for this plugin (must have less than " + Integer.MAX_VALUE + " points.");
			return;
		}
		
		String[] materialList = null;

		AmiraParameters parameters = null;
		if (AmiraParameters.isAmiraLabelfield(img)) {
			parameters = new AmiraParameters(img);
			materialList = parameters.getMaterialList();
		}

		ArrayList<Region> results = new ArrayList<Region>();

		ImageStack stack = img.getStack();

		byte[][] sliceDataBytes = null;
		float[][] sliceDataFloats = null;

		if (true) {
			sliceDataBytes = new byte[depth][];
			for (int z = 0; z < depth; ++z) {
				ByteProcessor bp = (ByteProcessor) stack.getProcessor(z+1);
				sliceDataBytes[z] = (byte[]) bp.getPixelsCopy();
			}
		} 
//		else {
//			sliceDataFloats = new float[depth][];
//			for (int z = 0; z < depth; ++z) {
//				FloatProcessor bp = (FloatProcessor) stack.getProcessor(z+1);
//				sliceDataFloats[z] = (float[]) bp.getPixelsCopy();
//			}
//		}
		
		// Preserve the calibration and color lookup tables
		// for generating new images of each individual
		// region.
//		Calibration calibration = img.getCalibration();

//		ColorModel cm = null;
//		if (ImagePlus.COLOR_256 == type) {
//			cm = stack.getColorModel();
//		}

//		ResultsTable rt=ResultsTable.getResultsTable();
//		rt.reset();

//		CancelDialog cancelDialog=new CancelDialog(this);
//		cancelDialog.show();

		boolean firstTime = true;

		int[]  labels 		= new int[depth * width * height];
		
		System.out.println("Before loop! There was "+results.size()+" regions.");
		
		while (true) {

//			if(pleaseStop)
//				break;
			
			/* Find one pixel that's above the minimum, or
			   find the maximum in the case where we're
			   not insisting that all regions are made up
			   of the same color.  These are set in all
			   cases... */
			
			int initial_x = -1;
			int initial_y = -1;
			int initial_z = -1;
			
			int foundValueInt = -1;
			float foundValueFloat = Float.MIN_VALUE;
			int maxValueInt = -1;
			float maxValueFloat = Float.MIN_VALUE;
			
			if (firstTime && startFromPointROI ) {
				
				initial_x = point_roi_x;
				initial_y = point_roi_y;
				initial_z = point_roi_z;

				if(byteImage)
					foundValueInt = sliceDataBytes[initial_z][initial_y * width + initial_x] & 0xFF;
				else
					foundValueFloat = sliceDataFloats[initial_z][initial_y * width + initial_x];

			} 
			
			if (true) { // byteImage && !startAtMaxValue

				// Just finding some point in the a region...
				for (int z = 0; z < depth && foundValueInt == -1; ++z) {
					for (int y = 0; y < height && foundValueInt == -1; ++y) {
						for (int x = 0; x < width; ++x) {
							int value = sliceDataBytes[z][y * width + x] & 0xFF;
							if (value > valuesOverDouble) {
								initial_x = x;
								initial_y = y;
								initial_z = z;
								foundValueInt = value;
								break;
							}
						}
					}
				}

				if (foundValueInt == -1) {
					break;
				}

			} 

			firstTime = false;

			int vint = foundValueInt;
			float vfloat = foundValueFloat;

			String materialName = null;
			if (materialList != null) {
				materialName = materialList[vint];
			}
			int pointsInQueue = 0;
			int queueArrayLength = 1024;
			int[] queue = new int[queueArrayLength];

			byte[] pointState 	= new byte[depth * width * height];
			
			
			int i = width * (initial_z * height + initial_y) + initial_x;
			pointState[i] = IN_QUEUE;
			queue[pointsInQueue++] = i;

			int pointsInThisRegion = 0;

			while (pointsInQueue > 0) {

				int nextIndex = queue[--pointsInQueue];

				int currentPointStateIndex = nextIndex;
				int pz = nextIndex / (width * height);
				int currentSliceIndex = nextIndex % (width * height);
				int py = currentSliceIndex / width;
				int px = currentSliceIndex % width;

				pointState[currentPointStateIndex] = ADDED;

//				if (byteImage) {
				sliceDataBytes[pz][currentSliceIndex] = 0;
//				} 
//				else {
//					sliceDataFloats[pz][currentSliceIndex] = Float.MIN_VALUE;
//				}
				
				++pointsInThisRegion;

				int x_unchecked_min = px - 1;
				int y_unchecked_min = py - 1;
				int z_unchecked_min = pz - 1;

				int x_unchecked_max = px + 1;
				int y_unchecked_max = py + 1;
				int z_unchecked_max = pz + 1;

				int x_min = (x_unchecked_min < 0) ? 0 : x_unchecked_min;
				int y_min = (y_unchecked_min < 0) ? 0 : y_unchecked_min;
				int z_min = (z_unchecked_min < 0) ? 0 : z_unchecked_min;

				int x_max = (x_unchecked_max >= width) ? width - 1 : x_unchecked_max;
				int y_max = (y_unchecked_max >= height) ? height - 1 : y_unchecked_max;
				int z_max = (z_unchecked_max >= depth) ? depth - 1 : z_unchecked_max;

				for (int z = z_min; z <= z_max; ++z) {
					for (int y = y_min; y <= y_max; ++y) {
						for (int x = x_min; x <= x_max; ++x) {

							// If we're not including diagonals,
							// skip those points.
							if ((!diagonal) && (x == x_unchecked_min || x == x_unchecked_max) && (y == y_unchecked_min || y == y_unchecked_max) && (z == z_unchecked_min || z == z_unchecked_max)) {
								continue;
							}
							int newSliceIndex = y * width + x;
							int newPointStateIndex = width * (z * height + y) + x;

							int neighbourValue = sliceDataBytes[z][newSliceIndex] & 0xFF;

							if (neighbourValue != vint) {
								continue;
							}

							if (0 == pointState[newPointStateIndex]) {
								pointState[newPointStateIndex] = IN_QUEUE;
								if (pointsInQueue == queueArrayLength) {
									int newArrayLength = (int) (queueArrayLength * 1.2);
									int[] newArray = new int[newArrayLength];
									System.arraycopy(queue, 0, newArray, 0, pointsInQueue);
									queue = newArray;
									queueArrayLength = newArrayLength;
								}
								queue[pointsInQueue++] = newPointStateIndex;
							}
						}
					}
				}
			}

			// So now pointState should have no IN_QUEUE
			// status points...
			Region region = new Region(vint, materialName, pointsInThisRegion, mustHaveSameValue);

			if (pointsInThisRegion < minimumPointsInRegionDouble) {
				// System.out.println("Too few points - only " + pointsInThisRegion);
				continue;
			}

			results.add(region);
			
			for (int z = 0; z < depth; ++z) {
				for (int y = 0; y < height; ++y) {
					for (int x = 0; x < width; ++x) {
						byte status = pointState[width * (z * height + y) + x];

						if (status == IN_QUEUE) {
							IJ.log("BUG: point " + x + "," + y + "," + z + " is still marked as IN_QUEUE");
						}

						if (status == ADDED) {
							//sliceBytes[y * width + x] = replacementValue;
							labels[width * (z * height + y) + x] = results.size();
						}
							
					}
				}
			}

			if ( (stopAfterNumberOfRegions > 0) && (results.size() >= stopAfterNumberOfRegions) ) {
				break;
			}
		}
		
		Collections.sort(results, Collections.reverseOrder());

		//System.out.println("Done! There was "+results.size()+" regions.");
		
		for (Iterator<Region> it = results.iterator(); it.hasNext();) {
			Region r = it.next();
			System.out.println(r.toString());		       			
		}
		
//		ImageStack newStack = new ImageStack(width, height);
//		ImagePlus out_labelled = NewImage.createRGBImage("Labelled_Connected_Regions", width, height, depth, NewImage.FILL_BLACK);
		// create table assigning each label a random color
		
		
		// this part has to be updated in the loop
		
		byte[][] label2color = randomColorSetGen(results.size());
		
		for (int z = 0; z < depth; ++z) {
			byte[] sliceBytesRed 	= new byte[width * height];
			byte[] sliceBytesGreen 	= new byte[width * height];
			byte[] sliceBytesBlue 	= new byte[width * height];
			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					
					sliceBytesRed[y * width + x] 	= label2color[][0]; // (byte)labels[width * (z * height + y) + x];
					sliceBytesGreen[y * width + x] 	= label2color[][1]; // (byte)labels[width * (z * height + y) + x];
					sliceBytesBlue[y * width + x] 	= label2color[][2]; // (byte)labels[width * (z * height + y) + x];
					
				}
			}
		}

//		ImagePlus newImagePlus = new ImagePlus("labels", newStack);
//		(new FileSaver(newImagePlus)).saveAsTiffStack("label.tif");
		
	}

	public void cancel() {
		pleaseStop = true;
	}
	
//	class CancelDialog extends Dialog implements ActionListener {
//		Button cancel;
//		Get_Connected_Regions plugin;
//		
//		public CancelDialog(Get_Connected_Regions plugin) {
//			super( IJ.getInstance(), "Find Connected Regions", false );
//			this.plugin = plugin;
//			cancel = new Button("Cancel 'Find Connected Regions'");
//			add(cancel);
//			cancel.addActionListener(this);
//			pack();
//		}
//	    
//		public void actionPerformed( ActionEvent e ) {
//			Object source = e.getSource();
//	                if( source == cancel ) {
//	                        plugin.cancel();
//			}
//		}
//		
//	}
	
	/* An inner class to make the results list sortable. */
 	private class Region implements Comparable {

		Region(int value, String materialName, int points, boolean sameValue) {
			byteImage = true;
			this.value = value;
			this.materialName = materialName;
			this.points = points;
			this.sameValue = sameValue;
		}

		Region(int points, boolean sameValue) {
			byteImage = false;
			this.points = points;
			this.sameValue = sameValue;
		}

		boolean byteImage;
		int points;
		String materialName;
		int value;
		boolean sameValue;

		public int compareTo(Object otherRegion) {
			Region o = (Region) otherRegion;
			return (points < o.points) ? -1 : ((points > o.points) ? 1 : 0);
		}

		@Override
		public String toString() {
			if (byteImage) {
				String materialBit = "";
				if (materialName != null) {
					materialBit = " (" + materialName + ")";
				}
				return "Region of value " + value + materialBit + " containing " + points + " points";
			} else {
				return "Region containing " + points + " points";
			}
		}

//		public void addRow( ResultsTable rt ) {
//			rt.incrementCounter();
//			if(byteImage) {
//				if(sameValue)
//					rt.addValue("Value in Region",value);
//				rt.addValue("Points In Region",points);
//				if(materialName!=null)
//					rt.addLabel("Material Name",materialName);
//			} else {
//				rt.addValue("Points in Region",points);
//			}
//		}

	}
 	
 	private byte[][] randomColorSetGen(int size){
 		// each row will be one random color
 		// random colors can be accessed with index
 		// same index will always give the same random color in rgb
 		byte[][] setCols = new byte[size][3];
 		
 		for (int i = 0; i < setCols.length; i++) {
 			setCols[i] = randomColorGen();
		}
 		return setCols;
 	}
 	
 	private byte[] randomColorGen(){
 		// returns r,g,b in bytes
 		Random generator = new Random();
 		
 		int randomIndex = generator.nextInt(256);
 		
 		byte[][] jet256 = new byte[][]{
 				{(byte)   0, (byte)   0, (byte) 131},  {(byte)   0, (byte)   0, (byte) 135},  
 				{(byte)   0, (byte)   0, (byte) 139},  {(byte)   0, (byte)   0, (byte) 143},  
 				{(byte)   0, (byte)   0, (byte) 147},  {(byte)   0, (byte)   0, (byte) 151},  
 				{(byte)   0, (byte)   0, (byte) 155},  {(byte)   0, (byte)   0, (byte) 159},  
 				{(byte)   0, (byte)   0, (byte) 163},  {(byte)   0, (byte)   0, (byte) 167},  
 				{(byte)   0, (byte)   0, (byte) 171},  {(byte)   0, (byte)   0, (byte) 175},  
 				{(byte)   0, (byte)   0, (byte) 179},  {(byte)   0, (byte)   0, (byte) 183},  
 				{(byte)   0, (byte)   0, (byte) 187},  {(byte)   0, (byte)   0, (byte) 191},  
 				{(byte)   0, (byte)   0, (byte) 195},  {(byte)   0, (byte)   0, (byte) 199},  
 				{(byte)   0, (byte)   0, (byte) 203},  {(byte)   0, (byte)   0, (byte) 207},  
 				{(byte)   0, (byte)   0, (byte) 211},  {(byte)   0, (byte)   0, (byte) 215},  
 				{(byte)   0, (byte)   0, (byte) 219},  {(byte)   0, (byte)   0, (byte) 223},  
 				{(byte)   0, (byte)   0, (byte) 227},  {(byte)   0, (byte)   0, (byte) 231},  
 				{(byte)   0, (byte)   0, (byte) 235},  {(byte)   0, (byte)   0, (byte) 239},  
 				{(byte)   0, (byte)   0, (byte) 243},  {(byte)   0, (byte)   0, (byte) 247},  
 				{(byte)   0, (byte)   0, (byte) 251},  {(byte)   0, (byte)   0, (byte) 255},  
 				{(byte)   0, (byte)   4, (byte) 255},  {(byte)   0, (byte)   8, (byte) 255},  
 				{(byte)   0, (byte)  12, (byte) 255},  {(byte)   0, (byte)  16, (byte) 255},  
 				{(byte)   0, (byte)  20, (byte) 255},  {(byte)   0, (byte)  24, (byte) 255},  
 				{(byte)   0, (byte)  28, (byte) 255},  {(byte)   0, (byte)  32, (byte) 255},  
 				{(byte)   0, (byte)  36, (byte) 255},  {(byte)   0, (byte)  40, (byte) 255},  
 				{(byte)   0, (byte)  44, (byte) 255},  {(byte)   0, (byte)  48, (byte) 255},  
 				{(byte)   0, (byte)  52, (byte) 255},  {(byte)   0, (byte)  56, (byte) 255},  
 				{(byte)   0, (byte)  60, (byte) 255},  {(byte)   0, (byte)  64, (byte) 255},  
 				{(byte)   0, (byte)  68, (byte) 255},  {(byte)   0, (byte)  72, (byte) 255},  
 				{(byte)   0, (byte)  76, (byte) 255},  {(byte)   0, (byte)  80, (byte) 255},  
 				{(byte)   0, (byte)  84, (byte) 255},  {(byte)   0, (byte)  88, (byte) 255},  
 				{(byte)   0, (byte)  92, (byte) 255},  {(byte)   0, (byte)  96, (byte) 255},  
 				{(byte)   0, (byte) 100, (byte) 255},  {(byte)   0, (byte) 104, (byte) 255},  
 				{(byte)   0, (byte) 108, (byte) 255},  {(byte)   0, (byte) 112, (byte) 255},  
 				{(byte)   0, (byte) 116, (byte) 255},  {(byte)   0, (byte) 120, (byte) 255},  
 				{(byte)   0, (byte) 124, (byte) 255},  {(byte)   0, (byte) 128, (byte) 255},  
 				{(byte)   0, (byte) 131, (byte) 255},  {(byte)   0, (byte) 135, (byte) 255},  
 				{(byte)   0, (byte) 139, (byte) 255},  {(byte)   0, (byte) 143, (byte) 255},  
 				{(byte)   0, (byte) 147, (byte) 255},  {(byte)   0, (byte) 151, (byte) 255},  
 				{(byte)   0, (byte) 155, (byte) 255},  {(byte)   0, (byte) 159, (byte) 255},  
 				{(byte)   0, (byte) 163, (byte) 255},  {(byte)   0, (byte) 167, (byte) 255},  
 				{(byte)   0, (byte) 171, (byte) 255},  {(byte)   0, (byte) 175, (byte) 255},  
 				{(byte)   0, (byte) 179, (byte) 255},  {(byte)   0, (byte) 183, (byte) 255},  
 				{(byte)   0, (byte) 187, (byte) 255},  {(byte)   0, (byte) 191, (byte) 255},  
 				{(byte)   0, (byte) 195, (byte) 255},  {(byte)   0, (byte) 199, (byte) 255},  
 				{(byte)   0, (byte) 203, (byte) 255},  {(byte)   0, (byte) 207, (byte) 255},  
 				{(byte)   0, (byte) 211, (byte) 255},  {(byte)   0, (byte) 215, (byte) 255},  
 				{(byte)   0, (byte) 219, (byte) 255},  {(byte)   0, (byte) 223, (byte) 255},  
 				{(byte)   0, (byte) 227, (byte) 255},  {(byte)   0, (byte) 231, (byte) 255},  
 				{(byte)   0, (byte) 235, (byte) 255},  {(byte)   0, (byte) 239, (byte) 255},  
 				{(byte)   0, (byte) 243, (byte) 255},  {(byte)   0, (byte) 247, (byte) 255},  
 				{(byte)   0, (byte) 251, (byte) 255},  {(byte)   0, (byte) 255, (byte) 255},  
 				{(byte)   4, (byte) 255, (byte) 251},  {(byte)   8, (byte) 255, (byte) 247},  
 				{(byte)  12, (byte) 255, (byte) 243},  {(byte)  16, (byte) 255, (byte) 239},  
 				{(byte)  20, (byte) 255, (byte) 235},  {(byte)  24, (byte) 255, (byte) 231},  
 				{(byte)  28, (byte) 255, (byte) 227},  {(byte)  32, (byte) 255, (byte) 223},  
 				{(byte)  36, (byte) 255, (byte) 219},  {(byte)  40, (byte) 255, (byte) 215},  
 				{(byte)  44, (byte) 255, (byte) 211},  {(byte)  48, (byte) 255, (byte) 207},  
 				{(byte)  52, (byte) 255, (byte) 203},  {(byte)  56, (byte) 255, (byte) 199},  
 				{(byte)  60, (byte) 255, (byte) 195},  {(byte)  64, (byte) 255, (byte) 191},  
 				{(byte)  68, (byte) 255, (byte) 187},  {(byte)  72, (byte) 255, (byte) 183},  
 				{(byte)  76, (byte) 255, (byte) 179},  {(byte)  80, (byte) 255, (byte) 175},  
 				{(byte)  84, (byte) 255, (byte) 171},  {(byte)  88, (byte) 255, (byte) 167},  
 				{(byte)  92, (byte) 255, (byte) 163},  {(byte)  96, (byte) 255, (byte) 159},  
 				{(byte) 100, (byte) 255, (byte) 155},  {(byte) 104, (byte) 255, (byte) 151},  
 				{(byte) 108, (byte) 255, (byte) 147},  {(byte) 112, (byte) 255, (byte) 143},  
 				{(byte) 116, (byte) 255, (byte) 139},  {(byte) 120, (byte) 255, (byte) 135},  
 				{(byte) 124, (byte) 255, (byte) 131},  {(byte) 128, (byte) 255, (byte) 128},  
 				{(byte) 131, (byte) 255, (byte) 124},  {(byte) 135, (byte) 255, (byte) 120},  
 				{(byte) 139, (byte) 255, (byte) 116},  {(byte) 143, (byte) 255, (byte) 112},  
 				{(byte) 147, (byte) 255, (byte) 108},  {(byte) 151, (byte) 255, (byte) 104},  
 				{(byte) 155, (byte) 255, (byte) 100},  {(byte) 159, (byte) 255, (byte)  96},  
 				{(byte) 163, (byte) 255, (byte)  92},  {(byte) 167, (byte) 255, (byte)  88},  
 				{(byte) 171, (byte) 255, (byte)  84},  {(byte) 175, (byte) 255, (byte)  80},  
 				{(byte) 179, (byte) 255, (byte)  76},  {(byte) 183, (byte) 255, (byte)  72},  
 				{(byte) 187, (byte) 255, (byte)  68},  {(byte) 191, (byte) 255, (byte)  64},  
 				{(byte) 195, (byte) 255, (byte)  60},  {(byte) 199, (byte) 255, (byte)  56},  
 				{(byte) 203, (byte) 255, (byte)  52},  {(byte) 207, (byte) 255, (byte)  48},  
 				{(byte) 211, (byte) 255, (byte)  44},  {(byte) 215, (byte) 255, (byte)  40},  
 				{(byte) 219, (byte) 255, (byte)  36},  {(byte) 223, (byte) 255, (byte)  32},  
 				{(byte) 227, (byte) 255, (byte)  28},  {(byte) 231, (byte) 255, (byte)  24},  
 				{(byte) 235, (byte) 255, (byte)  20},  {(byte) 239, (byte) 255, (byte)  16},  
 				{(byte) 243, (byte) 255, (byte)  12},  {(byte) 247, (byte) 255, (byte)   8},  
 				{(byte) 251, (byte) 255, (byte)   4},  {(byte) 255, (byte) 255, (byte)   0},  
 				{(byte) 255, (byte) 251, (byte)   0},  {(byte) 255, (byte) 247, (byte)   0},  
 				{(byte) 255, (byte) 243, (byte)   0},  {(byte) 255, (byte) 239, (byte)   0},  
 				{(byte) 255, (byte) 235, (byte)   0},  {(byte) 255, (byte) 231, (byte)   0},  
 				{(byte) 255, (byte) 227, (byte)   0},  {(byte) 255, (byte) 223, (byte)   0},  
 				{(byte) 255, (byte) 219, (byte)   0},  {(byte) 255, (byte) 215, (byte)   0},  
 				{(byte) 255, (byte) 211, (byte)   0},  {(byte) 255, (byte) 207, (byte)   0},  
 				{(byte) 255, (byte) 203, (byte)   0},  {(byte) 255, (byte) 199, (byte)   0},  
 				{(byte) 255, (byte) 195, (byte)   0},  {(byte) 255, (byte) 191, (byte)   0},  
 				{(byte) 255, (byte) 187, (byte)   0},  {(byte) 255, (byte) 183, (byte)   0},  
 				{(byte) 255, (byte) 179, (byte)   0},  {(byte) 255, (byte) 175, (byte)   0},  
 				{(byte) 255, (byte) 171, (byte)   0},  {(byte) 255, (byte) 167, (byte)   0},  
 				{(byte) 255, (byte) 163, (byte)   0},  {(byte) 255, (byte) 159, (byte)   0},  
 				{(byte) 255, (byte) 155, (byte)   0},  {(byte) 255, (byte) 151, (byte)   0},  
 				{(byte) 255, (byte) 147, (byte)   0},  {(byte) 255, (byte) 143, (byte)   0},  
 				{(byte) 255, (byte) 139, (byte)   0},  {(byte) 255, (byte) 135, (byte)   0},  
 				{(byte) 255, (byte) 131, (byte)   0},  {(byte) 255, (byte) 128, (byte)   0},  
 				{(byte) 255, (byte) 124, (byte)   0},  {(byte) 255, (byte) 120, (byte)   0},  
 				{(byte) 255, (byte) 116, (byte)   0},  {(byte) 255, (byte) 112, (byte)   0},  
 				{(byte) 255, (byte) 108, (byte)   0},  {(byte) 255, (byte) 104, (byte)   0},  
 				{(byte) 255, (byte) 100, (byte)   0},  {(byte) 255, (byte)  96, (byte)   0},  
 				{(byte) 255, (byte)  92, (byte)   0},  {(byte) 255, (byte)  88, (byte)   0},  
 				{(byte) 255, (byte)  84, (byte)   0},  {(byte) 255, (byte)  80, (byte)   0},  
 				{(byte) 255, (byte)  76, (byte)   0},  {(byte) 255, (byte)  72, (byte)   0},  
 				{(byte) 255, (byte)  68, (byte)   0},  {(byte) 255, (byte)  64, (byte)   0},  
 				{(byte) 255, (byte)  60, (byte)   0},  {(byte) 255, (byte)  56, (byte)   0},  
 				{(byte) 255, (byte)  52, (byte)   0},  {(byte) 255, (byte)  48, (byte)   0},  
 				{(byte) 255, (byte)  44, (byte)   0},  {(byte) 255, (byte)  40, (byte)   0},  
 				{(byte) 255, (byte)  36, (byte)   0},  {(byte) 255, (byte)  32, (byte)   0},  
 				{(byte) 255, (byte)  28, (byte)   0},  {(byte) 255, (byte)  24, (byte)   0},  
 				{(byte) 255, (byte)  20, (byte)   0},  {(byte) 255, (byte)  16, (byte)   0},  
 				{(byte) 255, (byte)  12, (byte)   0},  {(byte) 255, (byte)   8, (byte)   0},  
 				{(byte) 255, (byte)   4, (byte)   0},  {(byte) 255, (byte)   0, (byte)   0},  
 				{(byte) 251, (byte)   0, (byte)   0},  {(byte) 247, (byte)   0, (byte)   0},  
 				{(byte) 243, (byte)   0, (byte)   0},  {(byte) 239, (byte)   0, (byte)   0},  
 				{(byte) 235, (byte)   0, (byte)   0},  {(byte) 231, (byte)   0, (byte)   0},  
 				{(byte) 227, (byte)   0, (byte)   0},  {(byte) 223, (byte)   0, (byte)   0},  
 				{(byte) 219, (byte)   0, (byte)   0},  {(byte) 215, (byte)   0, (byte)   0},  
 				{(byte) 211, (byte)   0, (byte)   0},  {(byte) 207, (byte)   0, (byte)   0},  
 				{(byte) 203, (byte)   0, (byte)   0},  {(byte) 199, (byte)   0, (byte)   0},  
 				{(byte) 195, (byte)   0, (byte)   0},  {(byte) 191, (byte)   0, (byte)   0},  
 				{(byte) 187, (byte)   0, (byte)   0},  {(byte) 183, (byte)   0, (byte)   0},  
 				{(byte) 179, (byte)   0, (byte)   0},  {(byte) 175, (byte)   0, (byte)   0},  
 				{(byte) 171, (byte)   0, (byte)   0},  {(byte) 167, (byte)   0, (byte)   0},  
 				{(byte) 163, (byte)   0, (byte)   0},  {(byte) 159, (byte)   0, (byte)   0},  
 				{(byte) 155, (byte)   0, (byte)   0},  {(byte) 151, (byte)   0, (byte)   0},  
 				{(byte) 147, (byte)   0, (byte)   0},  {(byte) 143, (byte)   0, (byte)   0},  
 				{(byte) 139, (byte)   0, (byte)   0},  {(byte) 135, (byte)   0, (byte)   0},  
 				{(byte) 131, (byte)   0, (byte)   0},  {(byte) 128, (byte)   0, (byte)   0}   
 		};
 		
 		return jet256[randomIndex];
 	}
 	
 	private static final byte IN_QUEUE = 1;
	
 	private static final byte ADDED = 2;

}
