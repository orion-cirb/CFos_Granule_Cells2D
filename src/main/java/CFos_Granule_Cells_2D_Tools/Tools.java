package CFos_Granule_Cells_2D_Tools;

import CFos_Granule_Cells_2D_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import CFos_Granule_Cells_2D_Tools.Cellpose.CellposeTaskSettings;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.RGBStackMerge;
import ij.plugin.filter.Analyzer;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import static java.lang.Math.round;
import java.util.ArrayList;
import java.util.Collections;
import javax.swing.ImageIcon;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.Objects3DIntPopulationComputation;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.image3d.ImageHandler;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import org.apache.commons.io.FilenameUtils;


/**
 * @author ORION-CIRB
 */
public class Tools {
    private CLIJ2 clij2 = CLIJ2.getInstance();
    private final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    
    public String[] channelsName = {"Hoechst", "c-Fos"};
    public Calibration cal = new Calibration();
    public double pixArea = 0;
    
    // Nuclei detection
    private double meanNucArea = 50; 
    
    // c-Fos detection       
    private String cellposeEnvDirPath = (IJ.isWindows()) ? System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose" : "/opt/miniconda3/envs/cellpose";
    public String cellposeCFosModel = "cyto";
    private double resizeFactor = 0.5;
    public int cellposeCFosDiameter = 25;
    private double minCFosArea = 20;
    private double maxCFosArea = 80;
    private double cfosIntensityThresh = 250;
    private boolean useGpu = true;

    
    /**
     * Display a message in the ImageJ console and status bar
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }
    
    
    /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.showMessage("Error", "3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
    /**
     * Flush and close an image
     */
    public void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Find images extension
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        String[] files = imagesFolder.list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
            switch (fileExt) {
               case "nd" :
                   ext = fileExt;
                   break;
                case "czi" :
                   ext = fileExt;
                   break;
                case "lif"  :
                    ext = fileExt;
                    break;
                case "ics" :
                    ext = fileExt;
                    break;
                case "ics2" :
                    ext = fileExt;
                    break;
                case "lsm" :
                    ext = fileExt;
                    break;
                case "tif" :
                    ext = fileExt;
                    break;
                case "tiff" :
                    ext = fileExt;
                    break;
            }
        }
        return(ext);
    }
    
    
    /**
     * Find images in folder
     */
    public void findImages(String imagesFolder, String imageExt, ArrayList<String> imageFiles) {
        File inDir = new File(imagesFolder);
        File[] files = inDir.listFiles();
        
        for (File file: files) {
            if(file.isFile()) {
                String fileExt = FilenameUtils.getExtension(file.getName());
                if (fileExt.equals(imageExt) && !file.getName().startsWith("."))
                    imageFiles.add(file.getAbsolutePath());
            } else if (file.isDirectory() && !file.getName().equals("Results")) {
                findImages(file.getAbsolutePath(), imageExt, imageFiles);
            }
        }
        Collections.sort(imageFiles);
    }
    
    
    /**
     * Find image calibration
     */
    public void findImageCalib(IMetadata meta) {
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
    }
    
    
    /**
     * Find channels name
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws loci.common.services.DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                }
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                }
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelFluor(0, n);
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
                break;
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
                break;   
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        return(channels);         
    }
    
        
    /**
     * Generate dialog box
     */
    public String[] dialog(String[] channels) {      
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 160, 0);
        gd.addImage(icon);
        
        gd.addMessage("Channels", Font.getFont("Monospace"), Color.blue);
        int index = 0;
        for (String ch : channelsName) {
            gd.addChoice(ch, channels, channels[index]);
            index++;
        }

        gd.addMessage("Nuclei detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Mean nucleus area (µm2): ", meanNucArea);
        
        gd.addMessage("c-Fos cells detection", Font.getFont("Monospace"), Color.blue);
        gd.addDirectoryField("Cellpose environment path: ", cellposeEnvDirPath);
        
        gd.addNumericField("Min c-Fos area (µm2): ", minCFosArea);
        gd.addNumericField("Max c-Fos area (µm2): ", maxCFosArea);
        gd.addNumericField("Mean intensity threshold : ", cfosIntensityThresh);
        
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY pixel size (µm): ", cal.pixelWidth, 4);
        gd.showDialog();
        
        String[] chChoices = new String[channelsName.length];
        for (int n = 0; n < chChoices.length; n++) 
            chChoices[n] = gd.getNextChoice();
        if (gd.wasCanceled())
            chChoices = null;        
        
        meanNucArea = gd.getNextNumber();
        
        minCFosArea = gd.getNextNumber();
        maxCFosArea = gd.getNextNumber();
        cfosIntensityThresh = gd.getNextNumber();
        
        cal.pixelWidth = cal.pixelHeight = gd.getNextNumber();
        cal.pixelDepth = 1;
        pixArea = cal.pixelWidth*cal.pixelHeight;  
        
        return(chChoices);
    }
    
    
    /**
     * Get the number of nuclei in the image:
     * Compute nuclei total area and divide it by meanNucArea
     */
    public int getNbNuclei(ImagePlus img, Roi roi) {
        ImagePlus imgG = gaussian_filter(img, 4);
        imgG = threshold(imgG, "Otsu");
        imgG = median_filter(imgG, 4, 4);
        clearOutside(imgG, roi);
        imgG.setCalibration(cal);

        ResultsTable rt = new ResultsTable();
        Analyzer analyzer = new Analyzer(imgG, Analyzer.AREA+Analyzer.LIMIT, rt);
        IJ.setAutoThreshold(imgG, "Default dark");
        analyzer.measure();
        double area = rt.getValue("Area", 0);

        flush_close(imgG);
        return ((int) round(area/meanNucArea));
    }
     
    
    /**
     * 2D Gaussian filter using CLIJ2
     */ 
    public ImagePlus gaussian_filter(ImagePlus img, double sizeXY) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLGauss = clij2.create(imgCL);
       clij2.gaussianBlur2D(imgCL, imgCLGauss, sizeXY, sizeXY);
       clij2.release(imgCL);
       ImagePlus imgGauss = clij2.pull(imgCLGauss);
       clij2.release(imgCLGauss);
       return(imgGauss);
    }
    
    
    /**
     * Threshold using CLIJ2
     */
    public ImagePlus threshold(ImagePlus img, String thMed) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCLBin);
        return(imgBin);
    }
    
    
    /**
     * 2D median filter using CLIJ2
     */ 
    public ImagePlus median_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLMed = clij2.create(imgCL);
       clij2.median2DBox(imgCL, imgCLMed, sizeXY, sizeXY);
       clij2.release(imgCL);
       ImagePlus imgMed = clij2.pull(imgCLMed);
       clij2.release(imgCLMed);
       return(imgMed);
    }
    
    
    /**
     * Clear outside ROI
     */
    public void clearOutside(ImagePlus img, Roi roi) {
        PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), Roi.FREEROI);
        poly.setLocation(0, 0);
        for (int n = 1; n <= img.getNSlices(); n++) {
            ImageProcessor ip = img.getImageStack().getProcessor(n);
            ip.setRoi(poly);
            ip.setBackgroundValue(0);
            ip.setColor(0);
            ip.fillOutside(poly);
        }
        img.updateAndDraw();
    } 
    
    
    /*
     * Look for all cells in a 2D image with CellPose
     */
   public Objects3DIntPopulation cellposeDetection(ImagePlus img, Roi roi){
        // Resize image to speed up Cellpose computation
        ImagePlus imgResized = img.resize((int)(img.getWidth()*resizeFactor), (int)(img.getHeight()*resizeFactor), "none");

        // Define CellPose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(cellposeCFosModel, 1, cellposeCFosDiameter, cellposeEnvDirPath);
        settings.useGpu(useGpu);
        
        // Run CellPose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgResized);
        ImagePlus imgOut = cellpose.run(); 
        imgOut = imgOut.resize(img.getWidth(), img.getHeight(), "none");
        clearOutside(imgOut, roi);
        imgOut.setCalibration(cal);
        
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgOut));
        int nbCellsBeforeFiltering = pop.getNbObjects();
        System.out.println(nbCellsBeforeFiltering + " CellPose detections");
       
        
        
        Objects3DIntPopulationComputation popComputation = new Objects3DIntPopulationComputation​(pop);
        Objects3DIntPopulation popFilter = popComputation.getFilterSize(minCFosArea/pixArea, maxCFosArea/pixArea);
        filterDetectionsByIntensity(popFilter, img);
        popFilter.resetLabels();
        System.out.println(popFilter.getNbObjects() + " detections remaining after size and intensity filtering (" + (pop.getNbObjects()-popFilter.getNbObjects()) + " filtered out)");

        flush_close(imgResized);
        flush_close(imgOut);
        return(popFilter);
    }
   
   
    /**
     * Filter cells by intensity
     */
    public Objects3DIntPopulation filterDetectionsByIntensity(Objects3DIntPopulation cellPop, ImagePlus img) {
        cellPop.getObjects3DInt().removeIf(p -> 
                (new MeasureIntensity(p, ImageHandler.wrap(img)).getValueMeasurement(MeasureIntensity.INTENSITY_AVG) < cfosIntensityThresh));
        return(cellPop);
    }
    
    
    /**
     * Compute ROI area
     */
    public double roiArea(Roi roi, ImagePlus img) {
        PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), Roi.FREEROI);
        poly.setLocation(0, 0);
        img.setRoi(poly);
        img.setCalibration(cal);
        
        ResultsTable rt = new ResultsTable();
        Analyzer analyzer = new Analyzer(img, Analyzer.AREA, rt);
        analyzer.measure();
        return(rt.getValue("Area", 0));
    }   
    
    
    /*
     * Save population of cells in image
     */
    public void drawResults(Objects3DIntPopulation pop, ImagePlus img, String imageName, String outDir) {
        ImageHandler imgObj = ImageHandler.wrap(img).createSameDimensions();
        if (pop.getNbObjects() > 0)
            pop.drawInImage(imgObj);
        
        ImagePlus[] imgColors = {null, imgObj.getImagePlus(), null, img};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(cal);
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(outDir + imageName + ".tif"); 
        imgObj.closeImagePlus();
    }
    
}
