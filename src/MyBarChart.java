import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.general.DatasetUtilities;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/20/13
 * Time: 3:53 AM
 */
public class MyBarChart extends ApplicationFrame {

	public MyBarChart(String title, double[][] data) {
			super(title);

			final CategoryDataset dataset = createDataset(data);
			final JFreeChart chart = createChart(dataset);

			final ChartPanel chartPanel = new ChartPanel(chart);
			chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
			setContentPane(chartPanel);
	}

	private CategoryDataset createDataset(double[][] data) {

			return DatasetUtilities.createCategoryDataset("", "REGION", data);
	}

	private JFreeChart createChart(final CategoryDataset dataset) {

			final JFreeChart chart = ChartFactory.createBarChart("", "Category", "Region Average", dataset,
																		PlotOrientation.VERTICAL, true, true, false);
			return chart;
	}

//	public static void main(final String[] args) {
//					final double[][] data = new double[][] {
//														   {210,300,320,265,299},
//														   {200,304,201,201,340},
//			};
//		MyBarChart chart = new MyBarChart("Vertical Bar Chart Demo", data);
//			chart.pack();
//			RefineryUtilities.centerFrameOnScreen(chart);
//			chart.setVisible(true);
//	}

}
