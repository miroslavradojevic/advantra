package advantra.critpoint;

import java.awt.Button;
import java.awt.FlowLayout;
import java.awt.Panel;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Vector;

import javax.swing.JFileChooser;

import advantra.feature.MyHessian;
import advantra.processing.IntensityCalc;

import MLdetection.FeaturePool;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import imagescience.image.Axes;
import imagescience.image.ByteImage;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class TestPatches implements PlugIn, ActionListener {

	Button 		testImage, adaBoostConf;
	String 		test_file, conf_file;
	int 		patch_size;
	ImagePlus 	img;
	Image 		outimg;
	double 		s_start, s_end;
	int 		s_number;
	double[]	s;
	
	Image 		L1scales, L2scales, v1scales, v2scales;
	float[][]	v;
	
	public void run(String arg0) {
		
		test_file = Prefs.get("advantra.critpoint.TestPatches.test_file", "none");
		conf_file = Prefs.get("advantra.critpoint.TestPatches.conf_file", "none");
		
		s_start  = Prefs.get("advantra.critpoint.ExtractPatches.start_scale", 1.0);
		s_end    = Prefs.get("advantra.critpoint.ExtractPatches.end_scale", 8.0);
		s_number = Prefs.getInt("advantra.critpoint.ExtractPatches.nr_scales", 8);
		
		GenericDialog gd = new GenericDialog("Test Patches");
		
		Panel buttons_panel = new Panel();
		buttons_panel.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 1));
		
		testImage = new Button("Open test image");
		testImage.addActionListener(this);
		
		adaBoostConf = new Button("Open adaboost configuration");
		adaBoostConf.addActionListener(this);
		
		buttons_panel.add(testImage);
	    buttons_panel.add(adaBoostConf);
	    gd.addPanel(buttons_panel);
	    
	    gd.addMessage("default: "+test_file);
	    gd.addMessage("default: "+conf_file);
	    
	    gd.addMessage("---------------------\n");
	    
	    gd.addNumericField("start scale", s_start, 1);
		gd.addNumericField("end   scale", s_end, 1);
		gd.addNumericField("nr   scales", s_number, 0);
	    
		gd.showDialog();

	    if (gd.wasCanceled()) {
	        return;
	    }
	    
	    Prefs.set("advantra.critpoint.TestPatches.test_file", test_file);
	    Prefs.set("advantra.critpoint.TestPatches.conf_file", conf_file);
	    
	    s = new double[s_number];
		for (int i = 0; i < s_number; i++) {
			s[i] = (i==0)? s_start : s_start+i*((s_end-s_start)/(s_number-1));
		}
	    
	    // read conf file
	    File file = new File(conf_file);
	    if (!file.exists()) {
	        IJ.log("no adaboost config file");
	        return;
	    } 
	    
	    double[][] adaboost;
	    try {
	        BufferedReader br = new BufferedReader(new FileReader(file));
	        patch_size = Integer.valueOf(br.readLine());
	        int length = Integer.valueOf(br.readLine());
	        adaboost = new double[length][3];
	        for (int i = 0; i < length; i++) {
	            adaboost[i][0] = Double.valueOf(br.readLine());
	            adaboost[i][1] = Double.valueOf(br.readLine());
	            adaboost[i][2] = Double.valueOf(br.readLine());
	        }
	        br.close();
	    } catch (IOException e) {
	        return;
	    }
	    
	    IJ.log("adaboost parameters loaded");
	    
	    // create feature pool
	    FeaturePool featurePool = new FeaturePool(patch_size);
	    int fSize = featurePool.getSize();
	    IJ.log("feature pool size is " + fSize);
	    
	    // detection using the parameters
	    //go through the whole image and check all the patches
	    
	    
	    if(img!=null && (new File(conf_file)).exists()){
	    	img.setTitle("test image");
		    img.show();
		    IJ.open(conf_file);
	    }
	    else{
	    	return;
	    }
	    
	    Image inimg = Image.wrap(img);
	    Dimensions dims = inimg.dimensions();
	    outimg = new ByteImage(dims);
	    
	    v = new float[img.getWidth()*img.getHeight()][2];
	    
	    // extract hessian eigen values and eigen vectors for different scales
 		Dimensions new_dims 	= new Dimensions(img.getWidth(), img.getHeight(), s.length); // eigen analysis is stored per layer
	 		
	 	L1scales = new FloatImage(new_dims); L1scales.axes(Axes.X);
	 	L2scales = new FloatImage(new_dims); L2scales.axes(Axes.X);
	 	v1scales = new FloatImage(new_dims); v1scales.axes(Axes.X);
	 	v2scales = new FloatImage(new_dims); v2scales.axes(Axes.X);
	 		
	 	MyHessian my_hess = new MyHessian();
	 		
	 	double[] aL1 = new double[dims.x];
	 	double[] aL2 = new double[dims.x];
	 	double[] aV11 = new double[dims.x];
	 	double[] aV12 = new double[dims.x];
	 		
	 	Coordinates coords 	= new Coordinates();
	 		
	 	for (int i = 0; i < s.length; i++) {
	 		
	 		System.out.println("extracting hessian for scale "+IJ.d2s(s[i]));
	 		Vector<Image> hess = my_hess.eigs(inimg.duplicate(), s[i], true);
	 		// assign values to layers of L1, L2, v1, v2
	 		Image L2 	= hess.get(0); L2.axes(Axes.X);
	 		Image L1 	= hess.get(1); L1.axes(Axes.X);
	 		Image V11 	= hess.get(2); V11.axes(Axes.X);
	 		Image V12 	= hess.get(3); V12.axes(Axes.X);
	 			
	 		for (coords.y=0; coords.y<dims.y; ++coords.y) {
	 				
	 			coords.z = 0;
	 			L2.get(coords,aL2);
	 			L1.get(coords,aL1);
	 			V11.get(coords,aV11);
	 			V12.get(coords,aV12);
	 				
	 			coords.z = i;
	 			L2scales.set(coords, aL2);
	 			L1scales.set(coords, aL1);
	 			v1scales.set(coords, aV11);
	 			v2scales.set(coords, aV12);
	 				
	 		}
	 			
	 	}
	 	
	 	Coordinates cin = new Coordinates();
	    //Coordinates cout = new Coordinates();
	 	IntensityCalc im_calc = new IntensityCalc(img.getStack());
	 	//double[][] imageT;// = new double[patch_size][patch_size];
	 	int margin = (int)Math.ceil(patch_size/Math.sqrt(2));
	 	
	 	for (cin.x = margin; cin.x < dims.x-margin; cin.x++) {
			for (cin.y = margin; cin.y < dims.y-margin; cin.y++) {
				
				// extract direction (only at selected point)
				double[] vec = extractDirection(cin.x, cin.y);
				double theta = Math.atan(vec[1]/vec[0]);
//				float[][] locs = extractLocs(cin.x, cin.y, theta);
				
				double[][] imageT = extractPatch2D(cin.x, cin.y, theta, im_calc);
				
				int[][] imageIntT = getIntegralImage(imageT);
                int[] imFeaturesT = computeSelectedFeatures(imageIntT, featurePool, adaboost);
                int test = applyAdaBoost(adaboost, imFeaturesT);
                if (test == 1) {
//                    cout.x = cin.x + patch_size / 2;
//                    cout.y = cin.y + patch_size / 2;
//                    cout.t = cin.t;
                    outimg.set(cin, 255);
                } 
				
			}
		}
	 	
	 	
		IJ.log("directions extracted");
		outimg.imageplus().show();
	    // extract directions at each point inside margin
	        
	    IJ.log("done");
	    
	}

	public void actionPerformed(ActionEvent e) {
		if(e.getSource()==testImage){
			JFileChooser fc = new JFileChooser();
			int returnVal = fc.showOpenDialog(null);
	        if (returnVal != JFileChooser.APPROVE_OPTION) {
	            return;
	        }
	        
	        test_file = fc.getSelectedFile().getAbsolutePath();
	        IJ.showStatus("Opened " + test_file);
	        img = new ImagePlus(test_file);
	        
		}
		
		if(e.getSource()==adaBoostConf){
			JFileChooser fc = new JFileChooser();
			int returnVal = fc.showOpenDialog(null);
	        if (returnVal != JFileChooser.APPROVE_OPTION) {
	            return;
	        }
	        
	        conf_file = fc.getSelectedFile().getAbsolutePath();
	        IJ.showStatus("Opened " + conf_file);
	        
		}
	}

	private int[][] getIntegralImage(double[][] imageZ) {
	    int dimsx = imageZ[0].length;
	    int dimsy = imageZ.length;
	    int[][] imageInt = new int[dimsy][dimsx];
	    int[][] s = new int[dimsy][dimsx];
	    for (int j = 0; j < dimsy; j++) {
	        for (int i = 0; i < dimsx; i++) {
	            s[j][i] = (j - 1 >= 0) ? s[j - 1][i] + (int) imageZ[j][i] : (int) imageZ[j][i];
	            imageInt[j][i] = (i - 1 >= 0) ? imageInt[j][i - 1] + s[j][i] : s[j][i];
	        }
	    }
	    return imageInt;
	}
	
	private int[] computeSelectedFeatures(int[][] imageInt, FeaturePool featurePool, double[][] adaboost) {
//	    int dimsx = imageInt[0].length;
//	    int dimsy = imageInt.length;
	    int size = adaboost.length;
	    int[] imFeatures = new int[size];
	    for (int i = 0; i < size; i++) {
	        int j = (int) adaboost[i][0];
	        int[] r0 = featurePool.getFeature(j).getWhiteSquare();
	        int[] r1 = featurePool.getFeature(j).getBlackSquare();
	        double[] w = featurePool.getFeature(j).getWeights();

	        int i11 = imageInt[r1[1] + r1[3] - 1][r1[0] + r1[2] - 1];
	        int i00 = ((r1[0] - 1 < 0) || (r1[1] - 1 < 0)) ? 0 : imageInt[r1[1] - 1][r1[0] - 1];
	        int i01 = ((r1[0] - 1 < 0)) ? 0 : imageInt[r1[1] + r1[3] - 1][r1[0] - 1];
	        int i10 = ((r1[1] - 1 < 0)) ? 0 : imageInt[r1[1] - 1][r1[0] + r1[2] - 1];
	        int sum1 = i11 + i00 - i10 - i01;

	        i11 = imageInt[r0[1] + r0[3] - 1][r0[0] + r0[2] - 1];
	        i00 = ((r0[0] - 1 < 0) || (r0[1] - 1 < 0)) ? 0 : imageInt[r0[1] - 1][r0[0] - 1];
	        i01 = ((r0[0] - 1 < 0)) ? 0 : imageInt[r0[1] + r0[3] - 1][r0[0] - 1];
	        i10 = ((r0[1] - 1 < 0)) ? 0 : imageInt[r0[1] - 1][r0[0] + r0[2] - 1];
	        int sum0 = i11 + i00 - i10 - i01;

	        imFeatures[i] = (int) (w[0] * sum0 + w[1] * sum1);
	    }

	    return imFeatures;
	}
	
	int applyAdaBoost(double[][] adaboost, int[] imFeaturesT) {
		    double res = 0;
		    double object_res = 0;
		    for (int i = 0; i < adaboost.length; i++) {
		        object_res += adaboost[i][1];
		    }
		    object_res *= 0.5;

		    for (int i = 0; i < adaboost.length; i++) {
		        if (applyClassifier(imFeaturesT[i], adaboost[i][2])) {
		            res += adaboost[i][1];
		        }
		    }
		    int test = (res > object_res) ? 1 : 0;

		    return test;
	}
	 
	private boolean applyClassifier(double x, double thresh) {
		    return (x >= thresh) ? true : false;
	}
	
	private double[] 	extractDirection(int x, int y){

		double[] v = new double[2];
		
		Coordinates coord = new Coordinates();
		coord.x = x;	coord.y = y;
		
		L1scales.axes(Axes.Z);
		L2scales.axes(Axes.Z);
		
		double[] aL1 = new double[s.length];
		double[] aL2 = new double[s.length];
		
		L1scales.get(coord, aL1);
		L2scales.get(coord, aL2);
		
		double min_ratio = (aL2[0]<0)?Math.abs(aL1[0]/aL2[0]):Double.MAX_VALUE;
		int scale_idx = 0;
		
		for (int i = 0; i < s.length; i++) {
			if(aL2[i]<0){
				double ratio = Math.abs(aL1[i]/aL2[i]);
				if(ratio<min_ratio){
					min_ratio = ratio;
					scale_idx = i;
				}
			}
		}
		
		coord.z = scale_idx;
		v[0] = v1scales.get(coord);
		v[1] = v2scales.get(coord);
		
		return v;
	}
	
	private float[][] 	extractLocs(int x, int y, double theta){
		
		float[][] locs = new float[2][patch_size*patch_size]; 
		
		int N1 = patch_size/2;
		
		float x_beg = x-N1;
		float y_beg = y-N1;
		float x_cur, y_cur;
		
		int cnt = 0;
		
		for (int i = 0; i < patch_size; i++) {
			
			x_cur = x_beg + i;
			
			for (int j = 0; j < patch_size; j++) {
				
				y_cur = y_beg + j;
				
				locs[0][cnt] = x_cur-x;
				locs[1][cnt] = y_cur-y;
				cnt++;
				
			}
		}
		
		// rotate
		for (int i = 0; i < locs[0].length; i++) {
			float x_rot = locs[0][i]*(float)Math.cos(theta)-locs[1][i]*(float)Math.sin(theta);
			float y_rot = locs[0][i]*(float)Math.sin(theta)+locs[1][i]*(float)Math.cos(theta);
			locs[0][i] = x_rot;
			locs[1][i] = y_rot;
		}
		
		// return back
		for (int i = 0; i < locs[0].length; i++) {
			locs[0][i] += x;
			locs[1][i] += y;
		}
		
		return locs;
	}
	
	private byte[] 		extractPatch(int x, int y, double theta, IntensityCalc img_calc){
		
		byte[] 	patch 		= new byte		[patch_size*patch_size];
		float[][] 	locs 	= new float		[2][patch_size*patch_size];
		
		int N1 = patch_size/2;
		
		float x_beg = x-N1;
		float y_beg = y-N1;
		float x_cur, y_cur;
		
		int cnt = 0;
		for (int i = 0; i < patch_size; i++) {
			
			x_cur = x_beg + i;
			
			for (int j = 0; j < patch_size; j++) {
				
				y_cur = y_beg + j;
				
				locs[0][cnt] = x_cur-x;
				locs[1][cnt] = y_cur-y;
				cnt++;
				
			}
		}
		
		// rotate
		for (int i = 0; i < locs[0].length; i++) {
			float x_rot = locs[0][i]*(float)Math.cos(theta)-locs[1][i]*(float)Math.sin(theta);
			float y_rot = locs[0][i]*(float)Math.sin(theta)+locs[1][i]*(float)Math.cos(theta);
			locs[0][i] = x_rot;
			locs[1][i] = y_rot;
		}
		
		// return back
		for (int i = 0; i < locs[0].length; i++) {
			locs[0][i] += x;
			locs[1][i] += y;
		}
		
		cnt = 0;
		for (int i = 0; i < patch_size; i++) {
			for (int j = 0; j < patch_size; j++) {
				patch[cnt] = (byte)Math.round(img_calc.interpolateAt(locs[1][cnt], locs[0][cnt]));
				cnt++;
			}
		}
		
		return patch;
		
	}

	private double[][] 		extractPatch2D(int x, int y, double theta, IntensityCalc img_calc){
		
		double[][] 	patch 	= new double	[patch_size][patch_size];
		float[][] 	locs 	= new float		[2][patch_size*patch_size];
		
		int N1 = patch_size/2;
		
		float x_beg = x-N1;
		float y_beg = y-N1;
		float x_cur, y_cur;
		
		int cnt = 0;
		for (int i = 0; i < patch_size; i++) {
			
			x_cur = x_beg + i;
			
			for (int j = 0; j < patch_size; j++) {
				
				y_cur = y_beg + j;
				
				locs[0][cnt] = x_cur-x;
				locs[1][cnt] = y_cur-y;
				cnt++;
				
			}
		}
		
		// rotate
		for (int i = 0; i < locs[0].length; i++) {
			float x_rot = locs[0][i]*(float)Math.cos(theta)-locs[1][i]*(float)Math.sin(theta);
			float y_rot = locs[0][i]*(float)Math.sin(theta)+locs[1][i]*(float)Math.cos(theta);
			locs[0][i] = x_rot;
			locs[1][i] = y_rot;
		}
		
		// return back
		for (int i = 0; i < locs[0].length; i++) {
			locs[0][i] += x;
			locs[1][i] += y;
		}
		
		cnt = 0;
		for (int i = 0; i < patch_size; i++) {
			for (int j = 0; j < patch_size; j++) {
				patch[j][i] = Math.round(img_calc.interpolateAt(locs[1][cnt], locs[0][cnt])); // locs[0]~column, locs[1]~row
				//patch[cnt] = (byte)Math.round(img_calc.interpolateAt(locs[1][cnt], locs[0][cnt]));
				cnt++;
			}
		}
		
		return patch;
		
	}

	
}
