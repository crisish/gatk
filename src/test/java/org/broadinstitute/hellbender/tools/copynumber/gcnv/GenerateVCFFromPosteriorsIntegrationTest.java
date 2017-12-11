package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import org.broadinstitute.hellbender.CommandLineProgramTest;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.testng.annotations.Test;

/**
 * Created by asmirnov on 1/7/18.
 */
public class GenerateVCFFromPosteriorsIntegrationTest extends CommandLineProgramTest {

    // Test directory
    private static final File TEST_DIR = new File(toolsTestDir, "copynumber/gcnv");

    // gCNV calls output directory names
    private static final List<String> chunkDirectoriesNames = new ArrayList<>(Arrays.asList("/chunk_1/calls/", "/chunk_2/calls/"));
    private static final File REFERENCE_FILE = new File(hg19MiniReference);

    private static List<String> chunkDirectoriesFileList;

    static {
        chunkDirectoriesFileList = chunkDirectoriesNames.stream()
                .map(s -> TEST_DIR.getAbsolutePath().concat(s)).collect(Collectors.toList());
    }

    private static final String SAMPLE_NAME = "SAMPLE_0";

    @Test
    public void testGenerateVCF() {
        File outputVCF = createTempFile("test", ".vcf");

        final String chunkDirectoriesArguments = chunkDirectoriesFileList.stream()
                .reduce("-" + GenerateVCFFromPosteriors.POSTERIOR_CALL_DIRECTORY_FULL_NAME,
                        (s1, s2) -> s1 + " -" + GenerateVCFFromPosteriors.POSTERIOR_CALL_DIRECTORY_FULL_NAME + s2);

        final String[] arguments = {
                "-" + GenerateVCFFromPosteriors.POSTERIOR_CALL_DIRECTORY_FULL_NAME, chunkDirectoriesFileList.get(0),
                "-" + GenerateVCFFromPosteriors.POSTERIOR_CALL_DIRECTORY_FULL_NAME, chunkDirectoriesFileList.get(1),
                "-" + GenerateVCFFromPosteriors.SAMPLE_NAME_FULL_NAME, SAMPLE_NAME,
                "-" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, outputVCF.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE_FILE.getAbsolutePath()
        };

        runCommandLine(arguments);
        //do nothing
    }


}