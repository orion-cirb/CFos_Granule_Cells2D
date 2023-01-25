/*
 * In 2D images of granule cells, count nuclei and c-Fos-positive cells
 * Author: ORION-CIRB
 */

import CFos_Granule_Cells_2D_Tools.Tools;
import ij.*;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import mcib3d.geom2.Objects3DIntPopulation;
import org.apache.commons.io.FilenameUtils;


public class CFos_Granule_Cells_2D implements PlugIn {
    
    Tools tools = new Tools();
    private String imageDir = "";
    public String outDirResults = "";
    private BufferedWriter outPutResults;
    
    public void run(String arg) {
        try {
            if ((!tools.checkInstalledModules())) {
                return;
            }
            
            imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }
            
            // Find images with nd extension
            
            ArrayList<String> imageFiles = new ArrayList();
            tools.findImages(imageDir, "nd", imageFiles);
            if (imageFiles == null) {
                IJ.showMessage("Error", "No images found with nd extension");
                return;
            }
            
            // Determine if channels have tif or TIF extension
            String fileExt = tools.findImageType(new File(imageDir));
            if (imageFiles == null) {
                IJ.showMessage("Error", "No channels found with tif or TIF extension");
                return;
            }
            
            // Create output folder
            outDirResults = imageDir + File.separator + "Results" + File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            
            // Write header in results file
            String header = "Parent folder\tImage name\tROI name\tROI area (µm2)\tNb nuclei\tNb c-Fos cells\n";
            FileWriter  fwResults = new FileWriter(outDirResults + "results.xls", false);
            outPutResults = new BufferedWriter(fwResults);
            outPutResults.write(header);
            outPutResults.flush();
            
            // Find image calibration
            tools.cal.pixelWidth = tools.cal.pixelHeight = 0.1625;
            tools.cal.pixelDepth = 1;
            tools.pixArea = tools.cal.pixelWidth * tools.cal.pixelHeight * tools.cal.pixelDepth;
            
            // Find channels name
            String[] channels = {"w1CSU_405_t1" , "w4CSU_642_t1"};
            
            // Dialog box
            String[] chs = tools.dialog(channels);
            if (chs == null) {
                IJ.showMessage("Error", "Plugin canceled");
                return;
            }
            
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                String parentFolder = f.replace(imageDir, "").replace(FilenameUtils.getName(f), "");
                tools.print("--- ANALYZING IMAGE " + parentFolder + rootName + " ------");
                
                // Find ROI(s)
                String roiFile = imageDir+parentFolder+rootName+".roi";
                if (!new File(roiFile).exists())
                    roiFile = imageDir+parentFolder+rootName+".zip";
                if (!new File(roiFile).exists()) {
                    tools.print("ERROR: No ROI file found for image " + parentFolder + rootName);
                    IJ.showMessage("Error", "No ROI file found for image " + parentFolder + rootName);
                    continue;
                }
                
                RoiManager rm = new RoiManager(false);
                rm.runCommand("Open", roiFile);
                Roi[] rois = rm.getRoisAsArray();
                for (Roi roi : rois) {
                    String roiName = roi.getName();
                    tools.print("- Analyzing ROI " + roiName + " -");
                    
                    // Open Hoechst channel
                    tools.print("Opening nuclei channel...");
                    ImagePlus imgNuc = IJ.openImage(imageDir+parentFolder+rootName+"_"+chs[0]+"."+fileExt);
                    imgNuc.setRoi(roi);
                    ImagePlus imgNucCrop = imgNuc.crop();
                    tools.flush_close(imgNuc);
                    
                    // Compute nuclei number
                    tools.print("Counting nuclei...");
                    int nbNuclei = tools.getNbNuclei(imgNucCrop, roi);
                    System.out.println(nbNuclei + " nuclei found");
                    
                    // Open CFos channel
                    tools.print("Opening c-Fos cells channel...");
                    ImagePlus imgCFos = IJ.openImage(imageDir+parentFolder+rootName+"_"+chs[1]+"."+fileExt);
                    imgCFos.setRoi(roi);
                    ImagePlus imgCFosCrop = imgCFos.crop();
                    tools.flush_close(imgCFos);
                    
                    // Detect c-Fos cells with Cellpose
                    Objects3DIntPopulation cfosPop = tools.cellposeDetection(imgCFosCrop, roi);
                    int nbCFos = cfosPop.getNbObjects();
                    System.out.println(nbCFos+" c-Fos cells found");
                    
                    // Compute ROI area
                    double roiArea = tools.roiArea(roi, imgCFosCrop);
                    
                    // Write results
                    outPutResults.write(parentFolder.replace("/", "")+"\t"+rootName+"\t"+roiName+"\t"+roiArea+"\t"+nbNuclei+"\t"+nbCFos+"\n");
                    outPutResults.flush();
                    
                    // Save images
                    tools.drawResults(cfosPop, imgCFosCrop, parentFolder.replace("/", "_")+rootName+"_"+roiName, outDirResults);
                    
                    tools.flush_close(imgNucCrop);
                    tools.flush_close(imgCFosCrop);
                }
            }
            outPutResults.close();
            tools.print("--- All done! ---");
        } catch (IOException ex) {
            Logger.getLogger(CFos_Granule_Cells_2D.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
