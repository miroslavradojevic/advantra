package advantra.feature;

import mpicbg.imagefeatures.FloatArray2DSIFT;

public class MySIFT extends FloatArray2DSIFT {

	public MySIFT(final Param p){
		super(p);
	}
	
//	int width = img.getWidth();
//	int height = img.getHeight();
//	ImageProcessor ip1 = img.getProcessor().convertToFloat();
//	FloatArray2D fa = new FloatArray2D((float[])ip1.getPixels(), width, height);
//	FloatArray2DSIFT.Param p = new Param();
//	FloatArray2DSIFT sift = new FloatArray2DSIFT(p);	
//	sift.init(fa);
//	sift.run();
//	ArrayList<Feature> feat = new ArrayList<Feature>();
//	sift.extractFeatures(feat);
//	IJ.log("feat "+feat.size()+" features: ");
	
}
