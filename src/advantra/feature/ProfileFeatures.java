package advantra.feature;

import java.util.Vector;

public class ProfileFeatures {
	
	/*
	 * feature is a vector corresponding to directional response
	 * the aim of the feature is to count the number of peaks
	 * of the directional response, borders {bd} and filter levels
	 * example feature	: 
	 * ON {bd} OFF 
	 * ON{c1} <{b1}> OFF{c2} <{b2}> ON{c3} <{b3}> OFF{c4} 	: {} 
	 * 
	 */
	private Vector<Double> on;
	private Vector<Double> off;
	private Vector<Double> bd;
	
}
