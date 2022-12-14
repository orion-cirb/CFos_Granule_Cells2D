/*
 * In 2D images of granule cells, count nuclei
 * and count c-Fos-positive cells
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
    private boolean canceled = false;
    private String imageDir = "";
    public String outDirResults = "";
    private BufferedWriter outPutResults;
    
    
    public void run(String arg) {
        try {
            FileWriter fwResults = null;
            if (canceled) {
                IJ.showMessage("Plugin canceled");
                return;
            }
            if ((!tools.checkInstalledModules())) {
                return;
            }
            imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }
            // Find images with file_ext extension
            String file_ext = "nd";
            ArrayList<String> imageFiles = tools.findImages(imageDir, file_ext);
            if (imageFiles == null) {
                IJ.showMessage("Error", "No images found with " + file_ext + " extension");
                return;
            }
            // Create output folder
            outDirResults = imageDir + File.separator + "Results" + File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            // Write header in results file
            String header = "Image name\tROI name\tRoi area (µm2)\tNb nuclei\tNb c-Fos cells\n";
            fwResults = new FileWriter(outDirResults + "results.xls", false);
            outPutResults = new BufferedWriter(fwResults);
            outPutResults.write(header);
            outPutResults.flush();
            tools.cal.pixelWidth = tools.cal.pixelHeight = 0.1625;
            // dialog
            String[] chsName = {"405" , "642"};
            String[] channels = tools.dialog(chsName);
            if (channels == null) {
                IJ.showStatus("Plugin cancelled");
                return;
            }
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                // find rois
                String roiFile = imageDir+rootName+".roi";
                if (!new File(roiFile).exists())
                    roiFile = imageDir+rootName+".zip";
                if (!new File(roiFile).exists()) {
                    IJ.showMessage("Error", "No roi file found");
                    return;
                }
                
                // Find roi(s)
                RoiManager rm = new RoiManager(false);
                rm.runCommand("Open", roiFile);
                Roi[] rois = rm.getRoisAsArray();
                for (Roi roi : rois) {
                    String roiName = roi.getName();
                    tools.print("--- ANALYZING IMAGE " + rootName + " ------");
                    
                    // Open Hoechst channel
                    tools.print("- Analyzing DAPI channel -");
                    ImagePlus imgNuc = IJ.openImage(imageDir+rootName+"_w1CSU_405_t1.tif");
                    imgNuc.setRoi(roi);
                    ImagePlus imgNucCrop = imgNuc.crop();
                    tools.flush_close(imgNuc);
                    
                    // Find number of nuclei
                    System.out.println("Finding " + tools.channelNames[0] + " nuclei....");
                    double nbNuclei = tools.getNbNuclei(imgNucCrop, roi);
                    System.out.println(nbNuclei + " " + tools.channelNames[0] + " nuclei found in roi "+roiName);
                    
                    // Open CFos channel
                    tools.print("- Analyzing CFos channel -");
                    ImagePlus imgCFos = IJ.openImage(imageDir+rootName+"_w4CSU_642_t1.tif");
                    imgCFos.setRoi(roi);
                    ImagePlus imgCFosCrop = imgCFos.crop();
                    tools.flush_close(imgCFos);
                    Objects3DIntPopulation cfosPop = tools.cellPoseCellsPop(imgCFosCrop, roi);
                    System.out.println(cfosPop.getNbObjects()+" CFos cells found");
                    cfosPop = tools.filterDetectionsByIntensity(cfosPop, imgCFosCrop);
                    int nbCFos = cfosPop.getNbObjects();
                    System.out.println(nbCFos+" CFos cells found after threshold intensity in roi "+roiName);
                    
                    // Compute roi area
                    double roiArea = tools.roiArea(roi, imgCFosCrop);
                    
                    // Write results
                    outPutResults.write(rootName+"\t"+roiName+"\t"+roiArea+"\t"+nbNuclei+"\t"+nbCFos+"\n");
                    outPutResults.flush();
                    
                    // Save images
                    tools.saveImgObjects(cfosPop, imgNucCrop, imgCFosCrop, rootName+"_"+roiName, outDirResults);
                    tools.flush_close(imgNucCrop);
                    tools.flush_close(imgCFosCrop);
                }
            }
            outPutResults.close();
            System.out.println("Process done");
        } catch (IOException ex) {
            Logger.getLogger(CFos_Granule_Cells_2D.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
