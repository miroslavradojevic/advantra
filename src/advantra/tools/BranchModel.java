package advantra.tools;

import java.util.Random;

import flanagan.math.Matrix;

import advantra.general.ArrayHandling;
import advantra.general.DebugExport;
import advantra.general.ImageConversions;
import advantra.general.Transf;
import advantra.shapes.Sphere;
import advantra.shapes.Sphere.Planisphere_Extr_Mode;

import ij.ImagePlus;
import ij.io.FileSaver;

public class BranchModel { 		

	private static int 		MAX_RADIUS_STD_DEV 		= 16;
	private static int		MIN_HEIGHT 				= 32;
	private static int		MIN_WIDTH 				= 32;
	private static int		MIN_LENGTH 				= 8;
	
	// point coordinates p1 & p2 have fixed values for this 4-point model
	private double[] p1;
	private double[] p2;
	private double[] p3;
	private double[] p4;
	
	private byte[][] values;
	
	private int 		image_width;
	private int 		image_height;
	private int 		stack_size;
	private double 		radius_std_dev;
	
	public BranchModel(int height, int width, int length, double radius_std_dev){
		
		if(radius_std_dev>MAX_RADIUS_STD_DEV){
			System.err.println("Radius cannot be more than "+MAX_RADIUS_STD_DEV+" !");
			System.exit(1);
		}
		
		if(height<MIN_HEIGHT){
			System.err.println("Height cannot be less than "+MIN_HEIGHT+" !");
			System.exit(1);
		}
		
		if(width<MIN_WIDTH){
			System.err.println("Width cannot be less than "+MIN_WIDTH+" !");
			System.exit(1);
		}

		if(length<MIN_LENGTH){
			System.err.println("Length cannot be less than "+MIN_LENGTH+" !");
			System.exit(1);
		}
		
		/*
		 * allocate byte gray8 image to be used by methods
		 */
		
		
		this.image_width 	= width;
		this.image_height	= height;
		this.stack_size		= length;
		values = new byte[stack_size][image_width*image_height];
		this.radius_std_dev = radius_std_dev;

		System.out.println("Branch initialized with zero values... call drawModel() to complete the simulation.");
		
		
	} // constructor

	public void drawHorizontalModel(){



		p1 = new double[3];
		// fix the branch basis
		p1[0]		= image_height/8; 	// vertical position 		(row) 		(2d)
		p1[1]		= image_width/2; 	// horizontal position 		(column)	(2d)
		p1[2]		= stack_size/2; 	// layer
		p2 = new double[3];
		p2[0]		= image_height/2;
		p2[1]		= image_width/2;
		p2[2]		= stack_size/2;
		p3 = new double[3];
		p4 = new double[3];
		
		// side branches: 2 angles - chosen randomly - points p3, p4
		Random generator 		= new Random();
		double angle_branch_1 	= generator.nextDouble(); 		// uniformly distributed between 0 and 1
		angle_branch_1 			= angle_branch_1 * Math.PI/2; 	// rad wrt straight direction [striaght -> +pi/2]
		double angle_branch_2 	= generator.nextDouble();  
		angle_branch_2 			= angle_branch_2 * Math.PI/2;

		double branch_length	= Math.min(image_width/2, image_height/2)*0.7; // so that it fits inside
        System.out.print("here... ");
		createEndpoints(p1, p2, branch_length, angle_branch_1, angle_branch_2, p3, p4); // two endpoints at the end
		System.out.format("p1: %5.2f, %5.2f, %5.2f \n", p1[0], p1[1], p1[2]);
		System.out.format("p2: %5.2f, %5.2f, %5.2f \n", p2[0], p2[1], p2[2]);
		System.out.format("p3: %5.2f, %5.2f, %5.2f \n", p3[0], p3[1], p3[2]);
		System.out.format("p4: %5.2f, %5.2f, %5.2f \n", p4[0], p4[1], p4[2]);
		
		System.out.print("calling drawHorizontalModel()... ");
		
		values = new byte[stack_size][image_width*image_height];
		
		writeLineBetween(p1, p2);
		writeLineBetween(p2, p3);
		writeLineBetween(p2, p4);
		System.out.println("done.");
	}
	
	public void drawVerticalModel(){
		
		System.out.println("redefine points for vertical draw...");
		
		p1 = new double[3];
		p1[0]		= image_width/2; 	// vertical position 		(row) 		(2d)
		p1[1]		= image_height/2; 	// horizontal position 		(column)	(2d)
		p1[2]		= stack_size/8; 	// layer
		p2 = new double[3];
		p2[0]		= image_width/2;
		p2[1]		= image_height/2;
		p2[2]		= stack_size/2;
		p3 = new double[3];
		p4 = new double[3];
		
		// side branches: 2 angles - chosen randomly - points p3, p4
		Random generator 		= new Random();
		double angle_branch_1 	= generator.nextDouble(); 		// uniformly distributed between 0 and 1
		angle_branch_1 			= angle_branch_1 * Math.PI/2; 	// rad wrt straight direction [striaght -> +pi/2]
		double angle_branch_2 	= generator.nextDouble();  
		angle_branch_2 			= angle_branch_2 * Math.PI/2;

		double branch_length	= Math.min(Math.min((image_width/2)*0.7, (image_height/2)*0.7), (stack_size/2)*0.7); // so that it fits inside
		
		createEndpoints(p1, p2, branch_length, angle_branch_1, angle_branch_2, p3, p4); // two endpoints at the end

		System.out.format("p1: %5.2f, %5.2f, %5.2f \n", p1[0], p1[1], p1[2]);
		System.out.format("p2: %5.2f, %5.2f, %5.2f \n", p2[0], p2[1], p2[2]);
		System.out.format("p3: %5.2f, %5.2f, %5.2f \n", p3[0], p3[1], p3[2]);
		System.out.format("p4: %5.2f, %5.2f, %5.2f \n", p4[0], p4[1], p4[2]);
		
		
		// draw the lines between the points
		System.out.print("calling drawVerticalModel()... ");
		
		values = new byte[stack_size][image_width*image_height];
		
		writeLineBetween(p1, p2);
		writeLineBetween(p2, p3);
		writeLineBetween(p2, p4);
		
		System.out.println("done.");
	}
	
	public ImagePlus[] 	extractPlanispheres(double planisphere_radius, int Res, Planisphere_Extr_Mode mode){
		
		ImagePlus[] extracted_images = new ImagePlus[4];

		extracted_images[0] = extractPlanisphereAtPoint(p1, planisphere_radius, Res, mode); 
		extracted_images[1] = extractPlanisphereAtPoint(p2, planisphere_radius, Res, mode);
		extracted_images[2] = extractPlanisphereAtPoint(p3, planisphere_radius, Res, mode);
		extracted_images[3] = extractPlanisphereAtPoint(p4, planisphere_radius, Res, mode);
		
		return extracted_images;
		
	}
	
	public void 		exportPlanispheres( double planisphere_radius, int Res, Planisphere_Extr_Mode mode, String export_dir){ 
		
		// for p1-p4 points
		
		String[] file_names = new String[4];
		
		for (int i = 0; i < 4; i++) {
			file_names[i] = new String(export_dir+String.format("p%d_planisphere_r%2.1f_%dx%d_%s.tif", i, planisphere_radius, Res, 2*Res, mode));		
		}
		
		// p1
		(new FileSaver(extractPlanisphereAtPoint(p1, planisphere_radius, Res, mode))).saveAsTiff(file_names[0]); 
		System.out.format("saved to %s  \n", file_names[0]);
		// p2
		(new FileSaver(extractPlanisphereAtPoint(p2, planisphere_radius, Res, mode))).saveAsTiff(file_names[1]); 
		System.out.format("saved to %s  \n", file_names[1]);
		// p3
		(new FileSaver(extractPlanisphereAtPoint(p3, planisphere_radius, Res, mode))).saveAsTiff(file_names[2]); 
		System.out.format("saved to %s  \n", file_names[2]);
		// p4
		(new FileSaver(extractPlanisphereAtPoint(p4, planisphere_radius, Res, mode))).saveAsTiff(file_names[3]); 
		System.out.format("saved to %s  \n", file_names[3]);
		
	}
	
	public void 		exportPlanispheresAtPoints(double[][] point_rows, double planisphere_radius, int Res, Planisphere_Extr_Mode mode, String export_dir){
		
		int number_of_pts = point_rows.length;
		
		String[] names = new String[number_of_pts];
		
		for (int i = 0; i < number_of_pts; i++) {
			names[i] = new String(export_dir+String.format("x%.2f_y%.2f_z%.2f_r%.1f_%dx%d_%s.tif", 
					point_rows[i][0],
					point_rows[i][1],
					point_rows[i][2],
					planisphere_radius, 
					Res, 2*Res,
					mode
					));
			(new FileSaver(extractPlanisphereAtPoint(point_rows[i], planisphere_radius, Res, mode))).saveAsTiff(names[i]); 
			System.out.format("saved to %s  \n", names[i]);
		}
		
	}
	
	public void exportModelParams(String export_dir){
		String log_file_name 	= export_dir+"model_params.log";
		DebugExport metadata_file 	= new DebugExport(log_file_name);
		
		metadata_file.writeLine("width 				= "+image_width);
		metadata_file.writeLine("height 			= "+values[0].length/image_width);
		metadata_file.writeLine("length 			= "+values.length);
		metadata_file.writeLine("radius_std_dev		= "+radius_std_dev);
		
		double[] p;
		//p1
		p = getP1();	metadata_file.writeLine(String.format("p1: %2.2f, %2.2f, %2.2f", p[0], p[1], p[2]));
		//p2
		p = getP2();	metadata_file.writeLine(String.format("p2: %2.2f, %2.2f, %2.2f", p[0], p[1], p[2]));
		//p3
		p = getP3();	metadata_file.writeLine(String.format("p3: %2.2f, %2.2f, %2.2f", p[0], p[1], p[2]));
		//p4
		p = getP4();	metadata_file.writeLine(String.format("p4: %2.2f, %2.2f, %2.2f", p[0], p[1], p[2]));		
		
		System.out.println("log is recorded in:\t\t"+log_file_name);
		
		metadata_file.closeDebug();
		
	}
	
	private ImagePlus 	extractPlanisphereAtPoint(double[] planisphere_point, double planisphere_radius, int Res, Planisphere_Extr_Mode mode){
		
		Sphere sph = new Sphere(planisphere_point[0], planisphere_point[1], planisphere_point[2], planisphere_radius);
		
		return sph.extractPlanisphereView(ImageConversions.toGray8(values, image_width), Res, 0.7, mode);
		
	}
	
	private static void 	createEndpoints(
			double[] p_seed, 
			double[] p_branch, 
			double len, 
			double alfa_1, 
			double alfa_2, 
			double[] p_out_1,
			double[] p_out_2)
    {
		double[] Va = new double[3];
		Va[0] = p_branch[0] - p_seed[0];
		Va[1] = p_branch[1] - p_seed[1];
		Va[2] = p_branch[2] - p_seed[2];
		double V_norm = Math.sqrt(Va[0]*Va[0] + Va[1]*Va[1] + Va[2]*Va[2]);
		Va[0] = Va[0]/V_norm;
		Va[1] = Va[1]/V_norm;
		Va[2] = Va[2]/V_norm;
		
		double[] Vb = new double[3];
		double[] Vc = new double[3];
		
		Transf.cartesian(Va[0], Va[1], Va[2], Vb, Vc);

		Matrix V = new Matrix(
				new double[][]{
						{Va[0], Vb[0], Vc[0]},
						{Va[1], Vb[1], Vc[1]},
						{Va[2], Vb[2], Vc[2]}
						}
				);
		Matrix B  = Matrix.rowMatrix(
				new double[]{len*Math.cos(alfa_1), 0, len*Math.sin(alfa_1)}
				);
		
		Matrix A = B.times(Matrix.inverse(V));
		double[] a = A.getRowCopy(0);
		
		p_out_1[0] = a[0] + p_branch[0];
		p_out_1[1] = a[1] + p_branch[1];
		p_out_1[2] = a[2] + p_branch[2];
		
		B = Matrix.rowMatrix(
				new double[]{len*Math.cos(alfa_1), 0, -len*Math.sin(alfa_2)}
				);
		A = B.times(Matrix.inverse(V));
		a = A.getRowCopy(0);
		
		p_out_2[0] = a[0] + p_branch[0];
		p_out_2[1] = a[1] + p_branch[1];
		p_out_2[2] = a[2] + p_branch[2];
		
//		System.out.println("**** ****\n");
//		ArrayHandling.print1DArray(p_seed);
//		System.out.println("**** ****\n");
//		ArrayHandling.print1DArray(p_branch);
//		System.out.println("**** ****\n");
//		
//		ArrayHandling.print1DArray(new double[]{len*Math.cos(alfa_1), 0, len*Math.sin(alfa_1)});
//		System.out.println("**** ****\n");
//		ArrayHandling.print2DArray(new double[][]{
//						{Va[0], Vb[0], Vc[0]},
//						{Va[1], Vb[1], Vc[1]},
//						{Va[2], Vb[2], Vc[2]}
//						});
//		System.out.println("**** ****\n");
//		ArrayHandling.print1DArray(a);
	}
	
	private void writeLineBetween(
			double[] point_1, 
			double[] point_2
			//double std_dev_intensity,
			//byte[][] image_array, 
			//int image_width
			//boolean writeOnTopIfHigher
			){ 
		// static - same writing for every class
		
		int image_height 	= 	values[0].length/image_width;
		int image_len		=	values.length;
		
		// correct row (0-index) if they're out
		point_1[0] = (point_1[0]<0)?0:point_1[0];	point_1[0] = (point_1[0]>=image_height)?(image_height-1):point_1[0]; 
		point_2[0] = (point_2[0]<0)?0:point_2[0];	point_2[0] = (point_2[0]>=image_height)?(image_height-1):point_2[0];

		// correct col (1-index) if they're out
		point_1[1] = (point_1[1]<0)?0:point_1[1];	point_1[1] = (point_1[1]>=image_width )?(image_width -1):point_1[1];
		point_2[1] = (point_2[1]<0)?0:point_2[1];	point_2[1] = (point_2[1]>=image_width )?(image_width -1):point_2[1]; 
		
		// correct layer (2-index) if they're out
		point_1[2] = (point_1[2]<0)?0:point_1[2];	point_1[2] = (point_1[2]>=image_len)?(image_len-1):point_1[2]; 
		point_2[2] = (point_2[2]<0)?0:point_2[2];	point_2[2] = (point_2[2]>=image_len)?(image_len-1):point_2[2]; 
		
		/*
		 * unit orientation
		 */
		
		double[] n 	= new double[3];
		n[0] 		= point_2[0] - point_1[0];
		n[1] 		= point_2[1] - point_1[1];
		n[2] 		= point_2[2] - point_1[2];
		
		double 	n_len 	= Math.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
		
		if(n_len<=0.001){
			System.err.println("Error: distance p1-p2 was too small tocalculate direction: "+n_len);
			System.exit(1);
		}
		// normalize it
		n[0]        = n[0] / n_len;
		n[1]        = n[1] / n_len;
		n[2]        = n[2] / n_len;
		
		//go through every point
		for (int layer = 0; layer < values.length; layer++) {
			for (int coord = 0; coord < values[0].length; coord++) {
				
				int[] p = new int[3];
				p[0] = ArrayHandling.index2row_2d(coord, image_width);
				p[1] = ArrayHandling.index2col_2d(coord, image_width);
				p[2] = layer;

				/*
				 * check where within the tube it is
				 */
				
				// calculate (p1-p) dotprod (p2-p1)
				double projection_1 = 
						(point_1[0]-p[0])*(point_2[0]-point_1[0]) + 
						(point_1[1]-p[1])*(point_2[1]-point_1[1]) +
						(point_1[2]-p[2])*(point_2[2]-point_1[2]);
				
				// calculate p2-p
				double projection_2 = 
						(point_2[0]-p[0])*(point_1[0]-point_2[0]) +
						(point_2[1]-p[1])*(point_1[1]-point_2[1]) +
						(point_2[2]-p[2])*(point_1[2]-point_2[2]);
				
				double distance = 0;
				
				if(projection_1<=0 && projection_2<=0){

					//image_array[layer][ArrayHandling.sub2index_2d(p[0], p[1], image_width)] = (byte)255;
					// inside the line distance from line with a=p1 and n=unitVec(p1,p2)
					distance = Transf.distance_point_to_line_3d(point_1, n, p);
					
				}
				else if(projection_1>0 && projection_2<=0){
					
					// outside p1 - measure the distance to p1
					distance = dist(p, point_1);
					
				}
				else if(projection_1<=0 && projection_2>0){
					
					// outside p2 - measure the distance to p2
					distance = dist(p, point_2);
					
				}
				else{
					System.err.println("this case is not possible when printing line between p1 and p2");
					System.exit(1);
				}
				
				if(distance<3*radius_std_dev){//  // some reasonable distance from where we consider intensity is ~0
					
					byte value_to_write = (byte)Math.round(255 * Math.exp(-Math.pow(distance,2)/(2*Math.pow(radius_std_dev, 2))));
					
					values[layer][ArrayHandling.sub2index_2d(p[0], p[1], image_width)] = 
								(byte) Math.max(
										(int)(value_to_write & 0xff), 
										(int)(values[layer][ArrayHandling.sub2index_2d(p[0], p[1], image_width)] & 0xff)
										) ;
				}
			}
		}
	} // writeLineBetween

	private static double  dist(int[] p1, double[] p2){
		return Math.sqrt( Math.pow(p1[0]-p2[0], 2) + Math.pow(p1[1]-p2[1], 2) + Math.pow(p1[2]-p2[2], 2) );
	}
	
	public double[] getP1(){
		double[] p1_out = new double[3];
		p1_out[0]		= this.p1[0];
		p1_out[1]		= this.p1[1];
		p1_out[2]		= this.p1[2];
		return p1_out;
	}
	
	public double[] getP2(){
		double[] p2_out = new double[3];
		p2_out[0]		= this.p2[0];
		p2_out[1]		= this.p2[1];
		p2_out[2]		= this.p2[2];
		return p2_out;
	}
	
	public double[] getP3(){
		double[] p3_out = new double[3];
		p3_out[0]		= this.p3[0];
		p3_out[1]		= this.p3[1];
		p3_out[2]		= this.p3[2];
		return p3_out;
	}
	
	public double[] getP4(){
		double[] p4_out = new double[3];
		p4_out[0]		= this.p4[0];
		p4_out[1]		= this.p4[1];
		p4_out[2]		= this.p4[2];
		return p4_out;
	}

	public ImagePlus getModelAsImage(){
		return ImageConversions.toGray8(values, image_width);
	}
	
	public byte[][] getModelAsArray(){
		return this.values;
	}
	
	public void saveModel(String export_dir){
		
		(new FileSaver(ImageConversions.toGray8(values, image_width))).saveAsTiffStack(export_dir + String.format("branch_model.tif")); 
		System.out.println("exported : "+export_dir + String.format("branch_model.tif"));
		
	}
	
}