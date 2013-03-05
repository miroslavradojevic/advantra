package advantra.feature;

import java.io.File;
import java.io.IOException;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.CSVLoader;
import weka.core.converters.ConverterUtils.DataSink;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import imagescience.feature.Differentiator;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class Laplacian implements PlugIn {

	Image image;
	
	public void run(String arg0) {
		
		IJ.showMessage("Using Laplacian as a CP feature!");
		
		String image_path 	= System.getProperty("user.home")+File.separator+"train"+File.separator+"n01.tif";
		String csv_path		= System.getProperty("user.home")+File.separator+"train"+File.separator+"n01.csv";
		double sigma_start 	= 2.0;
		double sigma_end	= 3.0;
		int	sigma_nr 		= 2;
		
		GenericDialog gd = new GenericDialog("Testing laplacian...");
		gd.addStringField("image        :", image_path, 50);
		gd.addStringField("ground truth :", csv_path, 	50);
		gd.addNumericField("start 	sqrt(t) (sigma): ", sigma_start, 	2);
		gd.addNumericField("end 	sqrt(t) (sigma): ", sigma_end, 		2);
		gd.addNumericField("number of", 				sigma_nr, 		0);
		gd.showDialog();
		if (gd.wasCanceled()) return;
		image_path 	= gd.getNextString();
		sigma_start = (double)gd.getNextNumber();
		sigma_end 	= (double)gd.getNextNumber();
		sigma_nr 	= (int)gd.getNextNumber();
		
		ImagePlus 	img 	= new ImagePlus(image_path);
		Image 		image   = null;
		
		// calculate the features
		Differentiator df = new Differentiator();
		if(img!=null){
			if(img.getStackSize()>1){
				IJ.showMessage("Cannot work with 3d!");
				return;
			}
			image = new FloatImage(Image.wrap(img));
			//image.imageplus().show();
		}
		
		// scales
		double[] sigma = new double[sigma_nr];
		for (int i = 0; i < sigma_nr; i++) {
			sigma[i] = (i==0)? sigma_start : sigma_start+i*((sigma_end-sigma_start)/(sigma_nr-1));
		}
		
		// extract Lxx+Lyy at different scales
		Dimensions dims = image.dimensions();
		
		Image Lxx 		= new FloatImage(dims);
		Image Lyy		= new FloatImage(dims);
		Image Laplacian	= new FloatImage(dims);
		
		double[] aLxx 		= new double[dims.x];
		double[] aLyy 		= new double[dims.x];
		double[] aLaplacian = new double[dims.x];
		
		File csv_file = new File(csv_path);
		Instances data = loadFromCsvFile(csv_file);
		
		double[][] locations = new double[data.numInstances()][2];
		for (int i = 0; i < data.numInstances(); i++) {
			locations[i][0] = data.instance(i).value(0);
			locations[i][1] = data.instance(i).value(1);
		}
		data.setClassIndex(data.numAttributes()-1);
		
		// create dataset attributes
		FastVector attributes = new FastVector();
		
		for (int i = 0; i < sigma_nr; i++) {
			
			String label = String.format("laplacian_S%.2f", sigma[i]);
			attributes.addElement(new Attribute(label));
			
		}
		
		// add the class attribute
		attributes.addElement(data.classAttribute());
		
		Instances train = new Instances("train", attributes, locations.length);
		for (int i = 0; i < locations.length; i++) {
			train.add(new Instance(sigma_nr+1));
		}
		train.setClassIndex(train.numAttributes()-1);
		
		for (int i = 0; i < sigma.length; i++) { // run through scales
			
			System.out.println("processing scale sigma = "+sigma[i]+" ");
			Lxx 		= df.run(image.duplicate(), sigma[i], 2, 0, 0);
			Lyy 		= df.run(image.duplicate(), sigma[i], 0, 2, 0);
			
			Lxx.axes(Axes.X);
			Lyy.axes(Axes.X);
			Laplacian.axes(Axes.X);
			
			Coordinates coords 	= new Coordinates();
			for (coords.y=0; coords.y<dims.y; ++coords.y) {
				Lxx.get(coords,aLxx);
				Lyy.get(coords,aLyy);
				for (int x=0; x<dims.x; ++x){
					aLaplacian[x] = aLxx[x]+aLyy[x];
				}
				Laplacian.set(coords,aLaplacian);
			}
			
			System.out.println("adding values to train set");
			
			for (int k = 0; k < locations.length; k++) {
				
				int col = (int)Math.round(locations[k][0]);
				int row = (int)Math.round(locations[k][1]);
				
				// take the value from image with feats
				Coordinates coord = new Coordinates(row, col);
				double feat_value = Laplacian.get(coord);
				// store it in the train dataset
				train.instance(k).setValue(i, feat_value);
				train.instance(k).setValue(train.numAttributes()-1, data.instance(k).classValue());
				
			}
		}
		
		String train_path = System.getProperty("user.home")+File.separator+"train.csv";
		try {
			DataSink.write(train_path, train);
		} catch (Exception e) {
			e.printStackTrace();
		}
		IJ.showMessage("Train set saved to "+train_path);
		
	}

 	private Instances loadFromCsvFile(File f)  {
 		
 		CSVLoader loader 	= new CSVLoader();
 		Instances data;
 		
 		try {
 			loader.setSource(f);
			data = loader.getDataSet();
			return data;
		} 
 		catch (IOException e) {
			System.out.println("There was a problem loading "+f.getName()+" file. Returning null.");
			e.printStackTrace();
			return null;
		}

 	}
	
}
