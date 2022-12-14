package CFos_Granule_Cells_2D_Tools.Cellpose;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import static java.io.File.separatorChar;

import ij.IJ;
import CFos_Granule_Cells_2D_Tools.Cellpose.CellposeTaskSettings;


public class Cellpose {

    static void execute(List<String> options, CellposeTaskSettings settings, Consumer<InputStream> outputHandler) throws IOException, InterruptedException {
        List<String> cmd = new ArrayList<>();
        List<String> start_cmd = null ;

        // start terminal
        if (IJ.isWindows()) {
            start_cmd=  Arrays.asList("cmd.exe", "/C");
        } else if (IJ.isMacOSX() || IJ.isLinux()) {
            start_cmd = Arrays.asList("bash", "-c");
        }
        cmd.addAll(start_cmd);


        // Depending of the env type
        if (settings.envType.equals("conda")) {
            List<String> conda_activate_cmd = null;

            if (IJ.isWindows()) {
                // Activate the conda env
                conda_activate_cmd = Arrays.asList("CALL", "conda.bat", "activate", settings.cellposeEnvDirectory);
                cmd.addAll(conda_activate_cmd);
                // After starting the env we can now use cellpose
                cmd.add("&");// to have a second command
                List<String> cellpose_args_cmd = Arrays.asList("python", "-Xutf8", "-m", "cellpose");
                cmd.addAll(cellpose_args_cmd);
                // input options
                cmd.addAll(options);
            } else if (IJ.isMacOSX() || IJ.isLinux()) {
                // instead of conda activate (so much headache!!!) specify the python to use
                String python_path = settings.cellposeEnvDirectory+separatorChar+"bin"+separatorChar+"python";
                List<String> cellpose_args_cmd = new ArrayList<>(Arrays.asList( python_path , "-m","cellpose"));
                cellpose_args_cmd.addAll(options);

                // convert to a string
                cellpose_args_cmd = cellpose_args_cmd.stream().map(s -> {
                    if (s.trim().contains(" "))
                        return "\"" + s.trim() + "\"";
                    return s;
                }).collect(Collectors.toList());
                // The last part needs to be sent as a single string, otherwise it does not run
                String cmdString = cellpose_args_cmd.toString().replace(",","");

                // finally add to cmd
                cmd.add(cmdString.substring(1, cmdString.length()-1));
            }

        } else if (settings.envType.equals("venv")) { // venv

            if (IJ.isWindows()) {
                List<String> venv_activate_cmd = Arrays.asList("cmd.exe", "/C", new File(settings.cellposeEnvDirectory, "Scripts/activate").toString());
                cmd.addAll(venv_activate_cmd);
            } else if (IJ.isMacOSX() || IJ.isLinux()) {
                throw new UnsupportedOperationException("Mac/Unix not supported yet with virtual environment. Please try conda instead.");
            }

        } else {
            throw new UnsupportedOperationException("Virtual env type unrecognized!");
        }

        System.out.println(cmd.toString().replace(",", ""));
        ProcessBuilder pb = new ProcessBuilder(cmd).redirectErrorStream(true);

        Process p = pb.start();
        Thread t = new Thread(Thread.currentThread().getName() + "-" + p.hashCode()) {
            @Override
            public void run() {
                BufferedReader stdIn = new BufferedReader(new InputStreamReader(p.getInputStream()));
                try {
                    for (String line = stdIn.readLine(); line != null; ) {
                        System.out.println(line);
                        line = stdIn.readLine();// you don't want to remove or comment that line! no you don't :P
                    }
                } catch (IOException e) {
                    System.out.println(e.getMessage());
                }
            }
        };
        t.setDaemon(true);
        t.start();

        p.waitFor();

        int exitValue = p.exitValue();

        if (exitValue != 0) {
            System.out.println("Runner " + settings.cellposeEnvDirectory + " exited with value " + exitValue + ". Please check output above for indications of the problem.");
        } else {
            System.out.println(settings.envType + " , " + settings.cellposeEnvDirectory + " run finished");
        }

    }

}
