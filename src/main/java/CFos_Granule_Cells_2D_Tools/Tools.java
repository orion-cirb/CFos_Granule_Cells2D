package CFos_Granule_Cells_2D_Tools;

import CFos_Granule_Cells_2D_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import CFos_Granule_Cells_2D_Tools.Cellpose.CellposeTaskSettings;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import ij.plugin.filter.Analyzer;
import ij.process.ImageProcessor;
import io.scif.DependencyException;
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
import mcib3d.geom2.BoundingBox;
import mcib3d.geom2.Object3DInt;
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
    
    public String[] channelNames = {"Hoechst", "c-Fos"};
    public Calibration cal = new Calibration();
    public double pixArea = 0;
    
    // Nuclei detection
    private double meanNucArea = 50; 
    private double minCFosArea = 20;
    private double maxCFosArea = 100;
    
    // CFos detection       
    private String cellposeEnvDirPath = (IJ.isWindows()) ? System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose" : "/opt/miniconda3/envs/cellpose";
    public String cellposeCFosModel = "cyto";
    public int cellposeCFosDiameter = 25;
    private double cfosIntensityThresh = 200;
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
     * Find image type
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
                case "isc2" :
                    ext = fileExt;
                    break;
                case "tif" :
                    ext = fileExt;
                    break;    
            }
        }
        return(ext);
    }
    
    
    /**
     * Find images in folder
     */
    public ArrayList<String> findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in "+imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
    
    
    /**
     * Find image calibration
     */
    public Calibration findImageCalib(IMetadata meta) {
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
        return(cal);
    }
    
    
     /**
     * Find channels name
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
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
                        channels[n] = meta.getChannelName(0, n).toString();
                }
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelFluor(0, n).toString();
                break;
            case "ics" :
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
    public String[] dialog(String[] chs) {      
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 100, 0);
        gd.addImage(icon);
        gd.addMessage("Channels", Font.getFont("Monospace"), Color.blue);
        int index = 0;
        for (String chNames : channelNames) {
            gd.addChoice(chNames+" : ", chs, chs[index]);
            index++;
        }

        gd.addMessage("Nuclei detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Mean nucleus area (µm2): ", meanNucArea);
        
        gd.addMessage("c-Fos cells detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min CFos area (µm2): ", minCFosArea);
        gd.addNumericField("Max CFos area (µm2): ", maxCFosArea);
        gd.addNumericField("Intensity threshold : ", cfosIntensityThresh);
        gd.addDirectoryField("Cellpose environment path: ", cellposeEnvDirPath);
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY pixel size (µm): ", cal.pixelWidth);
        gd.showDialog();
        
        String[] chChoices = new String[channelNames.length];
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
     * Flush and close an image
     */
    public void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    /**
     * Clear out side roi
     * @param img
     * @param roi
     */
    public void clearOutSide(ImagePlus img, Roi roi) {
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
    
    /**
     * Get the number of nuclei in the image
     * compute nuclei area
     * divide area by meanNucleusArea
     */ 
    public double getNbNuclei(ImagePlus img, Roi roi) {
        ImagePlus imgG = gaussian_filter(img, 4);
        imgG = threshold(imgG, "Otsu");
        imgG = median_filter(imgG, 4, 4);
        clearOutSide(img, roi);
        IJ.setAutoThreshold(imgG, "Default dark");
        ResultsTable rt = new ResultsTable();
        Analyzer ana = new Analyzer(imgG, Analyzer.AREA+Analyzer.LIMIT, rt);
        ana.measure();
        double area = rt.getValue("Area", 0);
        rt.reset();
        img = imgG.duplicate();
        flush_close(imgG);
        return (round(area/meanNucArea));
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
     * Median filter using CLIJ2
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
     * Threshold 
     * USING CLIJ2
     * @param img
     * @param thMed
     */
    public ImagePlus threshold(ImagePlus img, String thMed) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCLBin);
        return(imgBin);
    }
    
    
    
    public Objects3DIntPopulation filterDetectionsByIntensity(Objects3DIntPopulation cellPop, ImagePlus img) {
        cellPop.getObjects3DInt().removeIf(p -> 
                (new MeasureIntensity(p, ImageHandler.wrap(img)).getValueMeasurement(MeasureIntensity.INTENSITY_AVG) < cfosIntensityThresh));
        cellPop.resetLabels();
        return(cellPop);
    }
    
    public double roiArea(Roi roi, ImagePlus img) {
        PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), Roi.FREEROI);
        poly.setLocation(0, 0);
        img.setRoi(poly);
        ResultsTable rt = new ResultsTable();
        Analyzer ana = new Analyzer(img, Analyzer.AREA, rt);
        ana.measure();
        double area = rt.getValue("Area", 0);
        return(area);
    }   
    
    
    /**
    * Find cells with cellpose
    * return cell cytoplasm
     * @param imgCell
     * @param roi
    * @return 
    */
    public Objects3DIntPopulation cellPoseCellsPop(ImagePlus imgCell, Roi roi){
        ImagePlus imgIn = null;
        clearOutSide(imgCell, roi);
        // resize to be in a friendly scale
        int width = imgCell.getWidth();
        int height = imgCell.getHeight();
        float factor = 0.5f;
        boolean resized = false;
        if (imgCell.getWidth() > 1024) {
            imgIn = imgCell.resize((int)(width*factor), (int)(height*factor), 1, "none");
            resized = true;
        }
        else
            imgIn = new Duplicator().run(imgCell);
        imgIn.setCalibration(cal);
        CellposeTaskSettings settings = new CellposeTaskSettings(cellposeCFosModel, 1, cellposeCFosDiameter, cellposeEnvDirPath);
        settings.useGpu(true);
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgIn);
        ImagePlus cellpose_img = cellpose.run(); 
        flush_close(imgIn);
        ImagePlus cells_img = (resized) ? cellpose_img.resize(width, height, 1, "none") : cellpose_img;
        flush_close(cellpose_img);
        cells_img.setCalibration(cal);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(cells_img));
        System.out.println(pop.getNbObjects());
        Objects3DIntPopulation cellFilterPop = new Objects3DIntPopulationComputation(pop).getFilterSize(minCFosArea/pixArea, maxCFosArea/pixArea);
        cellFilterPop.resetLabels();
        flush_close(cells_img);
        return(cellFilterPop);
    }

   
    
    /**
     * Label object
     * @param popObj
     * @param img 
     * @param fontSize 
     */
    public void labelObject(Object3DInt obj, ImagePlus img, int fontSize) {
        if (IJ.isMacOSX())
            fontSize *= 3;
        
        BoundingBox bb = obj.getBoundingBox();
        int z = bb.zmin + 1;
        int x = bb.xmin;
        int y = bb.ymin;
        img.setSlice(z);
        ImageProcessor ip = img.getProcessor();
        ip.setFont(new Font("SansSerif", Font.PLAIN, fontSize));
        ip.setColor(255);
        ip.drawString(Integer.toString((int)obj.getLabel()), x, y);
        img.updateAndDraw();
    }
        
     
     /**
     * Save foci Population in image
     * @param imgNuc nucleus blue channel
     * @param cFosPop CFos cells in green channel
     * @param imgCFos 
     * @param imageName 
     * @param outDir 
     */
    public void saveImgObjects(Objects3DIntPopulation cFosPop, ImagePlus imgNuc, ImagePlus imgCFos, String imageName, String outDir) {

        // Draw CFos cells in green
        ImageHandler imgObj = ImageHandler.wrap(imgCFos).createSameDimensions();
        if (cFosPop.getNbObjects() > 0)
            for (Object3DInt obj: cFosPop.getObjects3DInt())
                obj.drawObject(imgObj, obj.getLabel());
        
        // Save image
        ImagePlus[] imgColors = {null, imgObj.getImagePlus(), imgNuc};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(cal);
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(outDir + imageName + ".tif"); 
        imgObj.closeImagePlus();
    }
}
