package advantra.general;

import java.util.Vector;

public class Sort {

	public static void 		asc3   (double[] sort, int[] ind){
		double first   = sort[0]; 
		double second  = sort[1];
		double third   = sort[2]; 
		
		double a = Math.abs(first);
		double b = Math.abs(second);
		double c = Math.abs(third);
		
		if(a<b)
		{
			if(a<c)
			{
				if(b<c)
				{
					// a, b, c
					sort[0]=first; sort[1]=second; sort[2]=third;
					ind[0]=0; ind[1]=1; ind[2]=2;
				}
				else
				{
					// a, c, b
					sort[0]=first; sort[1]=third; sort[2]=second;
					ind[0]=0; ind[1]=2; ind[2]=1;
				}
			}
			else
			{
				// c, a, b
				sort[0]=third; sort[1]=first; sort[2]=second;
				ind[0]=2; ind[1]=0; ind[2]=1;				
			}
		}
		else
		{
			if(b<c)
			{
				if(a<c)
				{
					// b, a, c
					sort[0]=second; sort[1]=first; sort[2]=third;
					ind[0]=1; ind[1]=0; ind[2]=2;						
				}
				else
				{
					// b, c, a
					sort[0]=second; sort[1]=third; sort[2]=first;
					ind[0]=1; ind[1]=2; ind[2]=0;						
				}
			}
			else
			{
				// c, b, a
				sort[0]=third; sort[1]=second; sort[2]=first;
				ind[0]=2; ind[1]=1; ind[2]=0;				
			}
		}
	}
		
	public static double 	findMin(double[] in){
		double min_val = in[0];
		for (int i = 1; i < in.length; i++) {
			if(in[i]<min_val){
				min_val = in[i];
			}
		}
		return min_val;
	}
	
	public static float 	findMin(float[] in){
		float min_val = in[0];
		for (int i = 1; i < in.length; i++) {
			if(in[i]<min_val){
				min_val = in[i];
			}
		}
		return min_val;
	}
	
	public static double 	findMax(double[] in){
		double max_val = in[0];
		for (int i = 1; i < in.length; i++) {
			if(in[i]>max_val){
				max_val = in[i];
			}
		}
		return max_val;
	}	
	
	public static float 	findMax(float[] in){
		float max_val = in[0];
		for (int i = 1; i < in.length; i++) {
			if(in[i]>max_val){
				max_val = in[i];
			}
		}
		return max_val;
	}	
	
	public static float[] 	min_max(float[] a){
		
		float mn = a[0];
		float mx = a[0];
		
		for (int i = 1; i < a.length; i++) {
			if(a[i]<mn){
				mn = a[i];
			}
			if(a[i]>mx){
				mx = a[i];
			}
		}
		
		return new float[]{mn, mx};
	}

	public static double[] min_max(double[] a){
		
		double mn = a[0];
		double mx = a[0];
		
		for (int i = 1; i < a.length; i++) {
			if(a[i]<mn){
				mn = a[i];
			}
			if(a[i]>mx){
				mx = a[i];
			}
		}
		
		return new double[]{mn, mx};
	}
	
	public static double sum(double[] a){
		double out = a[0];
		for (int i = 1; i < a.length; i++) {
			out += a[i];
		}
		return out;
	}
	
	public static double[] sum_pos_neg(double[] a){
		double[] out = new double[2];
		for (int i = 0; i < a.length; i++) {
			if(a[i]>=0){
				out[0] += a[i];
			}
			else{
				out[1] += a[i];
			}
		}
		return out;
	}
	
	public static float sum(float[] a){
		float out = a[0];
		for (int i = 1; i < a.length; i++) {
			out += a[i];
		}
		return out;
	}
	
	public static Vector<Integer> localMax(float[] a, int N){
		
		Vector<Integer> out = new Vector<Integer>();
		
		for (int i = 0; i < a.length; i++) {
			
			float 	curr_max 		= Float.MIN_VALUE;
			int		curr_max_idx 	= 0;
			
			for (int j = i-(N/2); j <= i+(N/2); j++) {
				
				float read_val;
				if(j<0 || j>=a.length){
					read_val = Float.MIN_VALUE;
				}
				else{
					read_val = a[j];
				}
				
				if(read_val>curr_max){
					curr_max = read_val;
					curr_max_idx = j;
				}
				
			}
			
			if(i==curr_max_idx){
				out.add(i);
			}
			
		}
		
		return out;
	}
	
}
