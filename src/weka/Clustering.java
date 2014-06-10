package weka;

import weka.clusterers.SimpleKMeans;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;

import java.util.ArrayList;

/**
 * Created by miroslav on 6/9/14.
 */
public class Clustering {

	public static float[][] getKMeansDirectionsXY(ArrayList<float[]> _vxy, int _K) throws Exception {

		Instances directions2d = readDirections(_vxy);
		SimpleKMeans kMeans = kMeansDirections(directions2d, _K);

		Instances centroids = kMeans.getClusterCentroids();

		// get direction of every cluster
		float[][] out_directions = new float[_K][2];
		for (int i = 0; i < centroids.numInstances(); i++) {
			out_directions[i][0] = (float) centroids.instance(i).value(0);
			out_directions[i][1] = (float) centroids.instance(i).value(1);
		}

		return  out_directions;

	}

	private static Instances readDirections(ArrayList<float[]> _vxy) { // note: direction vxy have to be ||=1 unit vectors in order for the euclidean distance to work

		// declare attributes (vector coordinates)
		Attribute vx = new Attribute("vx");
		Attribute vy = new Attribute("vy");

		// declare feature vector (2d coordinate)
		FastVector att = new FastVector(2);
		att.addElement(vx);
		att.addElement(vy);

		// declare & fill train set up with instances
		Instances dirs_inst = new Instances("vxy", att, 10); // initialize with 10 instances
		for (int i = 0; i < _vxy.size(); i++) {
			// create the instance
			Instance iExample = new Instance(2);
			iExample.setValue((Attribute)att.elementAt(0), _vxy.get(i)[0]);
			iExample.setValue((Attribute)att.elementAt(1), _vxy.get(i)[1]);
			dirs_inst.add(iExample); // add the instance
		}

		return dirs_inst;
	}

	private static SimpleKMeans kMeansDirections(Instances _directions2d, int _K) {

		// create the model
		SimpleKMeans kMeans = new SimpleKMeans();
		try {
			kMeans.setNumClusters(_K);
			kMeans.buildClusterer(_directions2d);
		} catch (Exception e) {
			e.printStackTrace();
		}

		return kMeans;

	}

}
