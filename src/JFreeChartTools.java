import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.DefaultHighLowDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RectangleInsets;
import org.jfree.ui.RefineryUtilities;

import java.awt.*;
import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/20/13
 * Time: 3:28 PM
 */
public class JFreeChartTools extends ApplicationFrame {

    public JFreeChartTools(String title) {
        super(title);
    }

    // show XY coordinates
    public void JFreeChart_XYPlot(String title, String legend, ArrayList<double[]> Y)
    {
        //super(title);

        final XYSeries series = new XYSeries(legend);
        XYSeries seriesX = new XYSeries("x");
        XYSeries seriesY = new XYSeries("y");
        XYSeries seriesZ = new XYSeries("z");
        XYSeries seriesM = new XYSeries("m");

        for(double[] y : Y) {
            seriesX.add(y[3], y[0]);
            seriesY.add(y[3], y[1]);
            seriesZ.add(y[3], y[2]);
            seriesM.add(y[3], 1000*y[4]);
        }

        final XYSeriesCollection data = new XYSeriesCollection();
        data.addSeries(seriesX);
        data.addSeries(seriesY);
        data.addSeries(seriesZ);
        data.addSeries(seriesM);


        final JFreeChart chart = ChartFactory.createXYLineChart(
                title,
                "t",
                "Y",
                data,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        chart.setBackgroundPaint(Color.white);
        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));
        plot.setDomainGridlinePaint(Color.gray); //vertical grid lines
        plot.setRangeGridlinePaint(Color.gray); // horizontal grid lines

        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
        renderer.setSeriesPaint(3, Color.MAGENTA);
        renderer.setShapesVisible(true);
        renderer.setShapesFilled(true);

        //NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        //rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());

        final ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 370));
        setContentPane(chartPanel);
        pack();
        RefineryUtilities.centerFrameOnScreen(this);
        setVisible(true);
    }

    public void plotVecs(String title, String[] legends, double[][] vecs)
    {

        int nrPlots = legends.length;
        XYSeries[] series = new XYSeries[nrPlots];

        for (int p = 0; p < nrPlots; p++) {

            series[p] = new XYSeries(legends[p]);

            for (int i = 0; i < vecs[p].length; i++) {
                series[p].add(i, vecs[p][i]);
            }

        }

        final XYSeriesCollection data = new XYSeriesCollection(series[0]);
        for (int p = 1; p < nrPlots; p++) {
            data.addSeries(series[p]);
        }

        final JFreeChart chart = ChartFactory.createXYLineChart(
                title,
                "t",
                "Y",
                data,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        chart.setBackgroundPaint(Color.white);
        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));
        plot.setDomainGridlinePaint(Color.gray); //vertical grid lines
        plot.setRangeGridlinePaint(Color.gray); // horizontal grid lines

        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
        renderer.setSeriesShapesVisible(0, false);
        renderer.setSeriesPaint(0, Color.red);
        renderer.setSeriesShapesFilled(0, false);

        renderer.setSeriesShapesVisible(1, true);
        renderer.setSeriesPaint(1, Color.green);
        renderer.setSeriesShapesFilled(1, false);

        renderer.setSeriesShapesVisible(2, true);
        renderer.setSeriesPaint(2, Color.blue);
        renderer.setSeriesShapesFilled(2, true);


        final ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(1000, 370));
        setContentPane(chartPanel);
        pack();
        RefineryUtilities.centerFrameOnScreen(this);
        setVisible(true);
    }

    public void plotVec(String title, String xlabel, String ylabel, double[] Y)
    {
        final XYSeries series = new XYSeries(ylabel);

        for (int i = 0; i < Y.length; i++) {
            series.add(i, Y[i]);
        }

        final XYSeriesCollection data = new XYSeriesCollection(series);
        final JFreeChart chart = ChartFactory.createXYLineChart(
                title,
                xlabel,
                ylabel,
                data,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        chart.setBackgroundPaint(Color.white);
        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));
        plot.setDomainGridlinePaint(Color.gray); //vertical grid lines
        plot.setRangeGridlinePaint(Color.gray); // horizontal grid lines

        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
        renderer.setSeriesShapesVisible(1, true);
        renderer.setSeriesShapesFilled(1, false);
        final ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(1000, 370));
        setContentPane(chartPanel);
        pack();
        RefineryUtilities.centerFrameOnScreen(this);
        setVisible(true);
    }

    public void plotVec(String title, String xlabel, String ylabel, float[] Y)
    {
        final XYSeries series = new XYSeries(ylabel);

        double lowestLow = Double.MAX_VALUE;
        double highestHigh = Double.MIN_VALUE;

        for (int i = 0; i < Y.length; i++) {
            series.add(i, Y[i]);
            if (Y[i]<lowestLow)     lowestLow = Y[i];
            if (Y[i]>highestHigh)   highestHigh = Y[i];
        }

        final XYSeriesCollection data = new XYSeriesCollection(series);
        final JFreeChart chart = ChartFactory.createXYLineChart(
                title,
                xlabel,
                ylabel,
                data,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );


        chart.getXYPlot().getRangeAxis().setRange(lowestLow*0.95, highestHigh*1.05);

        chart.setBackgroundPaint(Color.white);
        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));
        plot.setDomainGridlinePaint(Color.gray); //vertical grid lines
        plot.setRangeGridlinePaint(Color.gray); // horizontal grid lines

        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
        renderer.setSeriesShapesVisible(1, true);
        renderer.setSeriesShapesFilled(1, false);
        final ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(600, 300));
        setContentPane(chartPanel);
        pack();
        RefineryUtilities.centerFrameOnScreen(this);
        setVisible(true);
    }


//    // plod displacements
//    public void JFCPlotDisplacements(String title, String legend, double[] Y)
//    {
//
//        final XYSeries series = new XYSeries(legend);
//
//        for (int i = 0; i < Y.length; i++) {
//            series.add(i, Y[i]);
//        }
//
//        final XYSeriesCollection data = new XYSeriesCollection(series);
//        final JFreeChart chart = ChartFactory.createXYLineChart(
//                title,
//                "t",
//                "Y",
//                data,
//                PlotOrientation.VERTICAL,
//                true,
//                true,
//                false
//        );
//
//        chart.setBackgroundPaint(Color.white);
//        XYPlot plot = (XYPlot) chart.getPlot();
//        plot.setBackgroundPaint(Color.white);
//        plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));
//        plot.setDomainGridlinePaint(Color.gray); //vertical grid lines
//        plot.setRangeGridlinePaint(Color.gray); // horizontal grid lines
//
//        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
//        renderer.setShapesVisible(true);
//        renderer.setShapesFilled(true);
//
//        //NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
//        //rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
//
//        final ChartPanel chartPanel = new ChartPanel(chart);
//        chartPanel.setPreferredSize(new java.awt.Dimension(800, 370));
//        setContentPane(chartPanel);
//        pack();
//        RefineryUtilities.centerFrameOnScreen(this);
//        setVisible(true);
//    }

/*    // plot histogram
    public void JFCPlotDisplacementHistogram(String title, String legend, double[] Y, double[] X)
    {
        // the next to lines to remove the gradient effect in displaying of bars
        ChartFactory.setChartTheme(StandardChartTheme.createLegacyTheme());
        BarRenderer.setDefaultBarPainter(new StandardBarPainter());

        HistogramDataset dataset = new HistogramDataset();
        dataset.setType(HistogramType.RELATIVE_FREQUENCY);
        dataset.addSeries("Histogram",Y, 20);
        dataset.addSeries("Histogram",X, 20);

        String plotTitle = "";
        String xaxis = "displacement";
        String yaxis = "Frequency";
        PlotOrientation orientation = PlotOrientation.VERTICAL;
        boolean show = false;
        boolean toolTips = false;
        boolean urls = false;
        JFreeChart chart = ChartFactory.createHistogram( plotTitle, xaxis, yaxis,
                dataset, orientation, show, toolTips, urls);

        chart.setBackgroundPaint(Color.white);

        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));
        plot.setDomainGridlinePaint(Color.gray); //vertical grid lines
        plot.setRangeGridlinePaint(Color.gray); // horizontal grid lines
        plot.setForegroundAlpha(0.6F);

        XYBarRenderer renderer = (XYBarRenderer) plot.getRenderer();
        renderer.setShadowVisible(false);
        renderer.setDrawBarOutline(true);
//        renderer.setShadowXOffset(0);
//        renderer.setShadowYOffset(0);
//        renderer.setShadowVisible(false);
//        renderer.setShapesFilled(true);

        //NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        //rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());

        final ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(450, 400));
        setContentPane(chartPanel);
        pack();
        RefineryUtilities.centerFrameOnScreen(this);
        setVisible(true);
    }

    public void plotHistogram(String title, String legend, ArrayList<Double> X, ArrayList<Double>Y, String path)
    {
        // the next to lines to remove the gradient effect in displaying of bars
        ChartFactory.setChartTheme(StandardChartTheme.createLegacyTheme());
        BarRenderer.setDefaultBarPainter(new StandardBarPainter());

        double[] xx = new double[X.size()];
        double[] yy = new double[Y.size()];
        for (int i = 0; i < X.size(); i++) {
            xx[i] = X.get(i);
        }
        for (int i = 0; i < Y.size(); i++) {
            yy[i] = Y.get(i);
        }

        HistogramDataset dataset = new HistogramDataset();
        dataset.setType(HistogramType.RELATIVE_FREQUENCY);
        //dataset.addSeries("Histogram",xx, 80); // immobile fraction
        dataset.addSeries("Histogram",yy, 80); // mobile fraction

        String plotTitle = title;
        String xaxis = "r";
        String yaxis = "Frequency";
        PlotOrientation orientation = PlotOrientation.VERTICAL;
        boolean show = false;
        boolean toolTips = false;
        boolean urls = false;
        JFreeChart chart = ChartFactory.createHistogram( plotTitle, xaxis, yaxis,
                dataset, orientation, show, toolTips, urls);

        chart.setBackgroundPaint(Color.white);

        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));
        plot.setDomainGridlinePaint(Color.gray); //vertical grid lines
        plot.setRangeGridlinePaint(Color.gray); // horizontal grid lines
        plot.setForegroundAlpha(0.5F);

        XYBarRenderer renderer = (XYBarRenderer) plot.getRenderer();
        renderer.setShadowVisible(false);
        renderer.setDrawBarOutline(true);

        renderer.setSeriesOutlinePaint(0, Color.white);
        //renderer.setSeriesPaint(0, Color.blue);
        renderer.setSeriesOutlinePaint(1, Color.white);
        //renderer.setSeriesPaint(1, Color.green);

        // add data from Marcel's histogram
        //load marcel's data
        final XYSeries series1 = new XYSeries("Random 1");
        try {
            BufferedReader br = new BufferedReader(new FileReader(
                    "/Users/stealth/data/Marcel/120703_ES-GG1_Cont_TypeII_TimeInc1_ParticleJumpHist.txt"));
            String str = br.readLine();
            String str1 = br.readLine();
            final StringTokenizer items = new StringTokenizer(str, "\t");
            final StringTokenizer items1 = new StringTokenizer(str1, "\t");

            for (int i = 0; i < 81; i++)
            {
                double r = Double.valueOf(items.nextToken()) * 1000;
                double x = Double.valueOf(items1.nextToken());
                series1.add(new Integer((int)r), new Double(x));
            }
            br.close();
        } catch (IOException ex) {
            ij.IJ.error("Unable to read particle positions (2)");
        }

        final XYSeriesCollection dataset2 = new XYSeriesCollection(series1);

        plot.setDataset(2, dataset2);
        XYAreaRenderer2 renderer2 = new XYAreaRenderer2();
        plot.setRenderer(2, renderer2);

        final ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(900, 700));
        setContentPane(chartPanel);
        pack();
        RefineryUtilities.centerFrameOnScreen(this);
        setVisible(true);
        try {
            ChartUtilities.saveChartAsJPEG(new File(path), chart, 900, 700);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }*/

}
