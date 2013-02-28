package advantra.critpoint;

import java.io.File;

import weka.core.Instances;
import weka.core.Utils;
import weka.core.converters.CSVLoader;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Remove;

public class ExtractFeatures   { // implements PlugIn

//	public void run(String arg0) {
//		
//		
//	}
	
	public static void main(String[] args) throws Exception {
		
		File f 				= new File(args[0]);
		double[][] f_data 	= extractCols(f, new int[]{0, 1});
		System.out.println(Utils.arrayToString(f_data));
		
//		CSVLoader loader = new CSVLoader();
//		loader.setSource(new File(args[0]));
//		Instances gnd_tth = loader.getDataSet();
//		
//		// remove last attribute 
//		Remove rm = new Remove();
//		rm.setAttributeIndices(""+gnd_tth.numAttributes());
//		rm.setInputFormat(gnd_tth);
//		Instances gnd_tth_locs = Filter.useFilter(gnd_tth, rm);
//		
//		System.out.println(gnd_tth_locs);
//		
//		double[][] a = new double[gnd_tth_locs.numInstances()][gnd_tth_locs.numAttributes()];
//		
//		for (int i = 0; i < gnd_tth_locs.numInstances(); i++) {
//			a[i] = gnd_tth_locs.instance(i).toDoubleArray();
//		}
		
		
		
//		Instances data0 = DataSource.read( (new File(args[0])).getAbsolutePath() );
//		Instances data1 = DataSource.read( (new File(args[1])).getAbsolutePath() );
//		
//		System.out.println("data0:\n"+data0.numInstances()+" x "+data0.numAttributes());
//		System.out.println("data1:\n"+data1.numInstances()+" x "+data1.numAttributes());
//		
//		String[] options = Utils.splitOptions("help");//new String[]{"append"};
//		
//		Instances.main(options);
		
		//Instances data  = Instances.mergeInstances(data0, data1);
		
		//System.out.println(" *** AFTER MERGING: *** \n\n"+data);
	}
	
	private static double[][] extractCols(File csv_file, int[] attribs_to_keep) throws Exception {
		
		CSVLoader loader = new CSVLoader();
		loader.setSource(csv_file);
		Instances data = loader.getDataSet();
		
		// remove last attribute 
		Remove rm = new Remove();
		rm.setInvertSelection(true); // columns are kept
		rm.setAttributeIndicesArray(attribs_to_keep);//""+data.numAttributes()+"");
		rm.setInputFormat(data);
		Instances data_no_last_att = Filter.useFilter(data, rm);
		
		// convert to double[][]
		
		double[][] out = new double[data_no_last_att.numInstances()][data_no_last_att.numAttributes()];
		
		for (int i = 0; i < data_no_last_att.numInstances(); i++) {
			out[i] = data_no_last_att.instance(i).toDoubleArray();
		}
		
		return out;
		
	}

}
