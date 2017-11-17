package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.commons.io.output.NullOutputStream;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.ReadLengthReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.utils.ClassUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.mockito.internal.util.collections.Sets;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.PrintStream;
import java.util.*;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class GATKAnnotationPluginDescriptorUnitTest {
    // null print stream for the tests
    private static final PrintStream nullMessageStream = new PrintStream(new NullOutputStream());

//======================================================================================================================
//Methods from InbreedingCoefficient for testing purposes
    private static final Allele Aref = Allele.create("A", true);
    private static final Allele T = Allele.create("T");
    private static final Allele C = Allele.create("C");

    //make sure that compound hets (with no ref) don't add to het count
    VariantContext inbreedingCoefficientVC = makeVC("1", Arrays.asList(Aref, T, C),
            makeG("s1", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s2", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s3", T, C, 7099, 2530, 7099, 3056, 0, 14931),
            makeG("s4", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s5", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s6", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s7", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s8", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s9", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s10", Aref, T, 2530, 0, 7099, 366, 3056, 14931),

            //add a bunch of hom samples that will be ignored if we use s1..s10 as founders
            makeG("s11", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s12", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s13", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s14", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s15", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s16", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s17", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931)
    );
    private VariantContext makeVC(String source, List<Allele> alleles, Genotype... genotypes) {
        int start = 10;
        int stop = start; // alleles.contains(ATC) ? start + 3 : start;
        return new VariantContextBuilder(source, "1", start, stop, alleles).genotypes(Arrays.asList(genotypes)).filters((Set<String>) null).make();
    }
    private Genotype makeG(String sample, Allele a1, Allele a2, int... pls) {
        return new GenotypeBuilder(sample, Arrays.asList(a1, a2)).PL(pls).make();
    }
//======================================================================================================================

    private List<Annotation> instantiateFilter(final CommandLineParser clp) {
        GATKAnnotationPluginDescriptor readFilterPlugin = clp.getPluginDescriptor(GATKAnnotationPluginDescriptor.class);
        return Arrays.asList(readFilterPlugin.getMergedAnnotations().toArray(new Annotation[0]));
    }

    @DataProvider
    public Object[][] badAnnotationGroupsDataProvider() {
        Object[][] out = {
                { Arrays.asList("StandardAnnotation", "Foo", "Bar")},
                { Arrays.asList("RMSMappingQuality")},
                { Arrays.asList("Annotation")},
                { Arrays.asList("Object")}};
        return out;
    }

    @Test (dataProvider = "badAnnotationGroupsDataProvider", expectedExceptions = GATKException.class)
    public void testInvalidRequestedAnnotationGroup(List<String> testGroups) { //TODO this will probably turn into a check for class inheretance
        //This test asserts that the plugin descriptior will crash if an invalid annotation group is requested
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(null, testGroups)),
                Collections.emptySet());
        String[] args = {};
        clp.parseArguments(nullMessageStream, args);
        Collection<Annotation> annots = instantiateFilter(clp);
    }

    @DataProvider
    public Object[][] badAnnotationsDataProvider() {
        Object[][] out = {
                { Arrays.asList("-A", "StandardAnnotation")},
                { Arrays.asList("-A", "RMSMappingQuality", "-A", "RMSMappingQuality")},
                { Arrays.asList("-A", "RMSMappingQuality", "-AX", "RMSMappingQuality")},
                { Arrays.asList("-A", "foo")},
                { Arrays.asList("-A", "VariantAnnotator")}};
        return out;
    }

    @Test (dataProvider = "badAnnotationsDataProvider", expectedExceptions = CommandLineException.class)
    public void testInvalidRequestedAnnotations(List<String> arguments) { //TODO this will probably turn into a check for class inheretance
        //This test asserts that the plugin descriptior will crash if an invalid annotation group is requested
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(null, null)),
                Collections.emptySet());
        String[] args = arguments.toArray(new String[arguments.size()]);
        clp.parseArguments(nullMessageStream, args);
        Collection<Annotation> annots = instantiateFilter(clp);
    }

    @DataProvider
    public Object[][] annotationsWithRequiredArguments(){
        return new Object[][]{{ InbreedingCoeff.class.getSimpleName(), "--founderID", "s1"}};
    }

    // fail if a filter with required arguments is specified without corresponding arguments
    // TODO this test is disabled as there are currently no annotations for which the argument is required
    @Test(dataProvider = "annotationsWithRequiredArguments", expectedExceptions = CommandLineException.MissingArgument.class, enabled = false)
    public void testDependentFilterArguments(
            final String annot,
            final String argName,   //unused
            final String[] argValue) { //unused
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(null, null)),
                Collections.emptySet());

        String[] args = {"--annotation", annot};  // no args, just enable filters
        clp.parseArguments(nullMessageStream, args);
    }

    // fail if a filter's arguments are passed but the filter itself is not enabled
    @Test(dataProvider = "annotationsWithRequiredArguments", expectedExceptions = CommandLineException.class)
    public void testDanglingAnnotationArguments(
            final String filter,
            final String argName,
            final String argValue) {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(null, null)),
                Collections.emptySet());

        String[] args = { argName, argValue }; // no read filter set

        // no need to instantiate the filters - dependency errors are caught by the command line parser
        clp.parseArguments(nullMessageStream, args);
    }

    //Filters with arguments to verify filter test method actually filters
    @DataProvider
    public Object[][] annotationsWithGoodArguments(){
        return new Object[][]{
                { InbreedingCoeff.class.getSimpleName(), (Consumer<Annotation>)a -> {
                            Assert.assertEquals(Double.valueOf((String) ((InbreedingCoeff) a)
                                            .annotate(null, inbreedingCoefficientVC, null)
                                            .get(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY)),
                                    -0.3333333, 0.001, "InbreedingCoefficientScores"); },
                        "--founderID", new String[]{"s1", "s2",  "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10"} },};
    }

    @Test(dataProvider = "annotationsWithGoodArguments")
    public void testAnnotationArguments(
            final String annotation,
            final Consumer<Annotation> condition,
            final String argName,
            final String[] argValues) {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(null, null)),
                Collections.emptySet());

        List<String> args = Stream.of("--annotation", annotation).collect(Collectors.toList());
        Arrays.asList(argValues).forEach(arg -> {args.addAll(Arrays.asList(argName, arg));});

        clp.parseArguments(nullMessageStream, args.toArray(new String[args.size()]));
        List<Annotation> rf = instantiateFilter(clp);
        Assert.assertEquals(rf.size(), 1);
        condition.accept(rf.get(0));
    }


    @Test
    public void testToolDefaultAnnotationArgumentsOverriding() { //TODO this is a hot mess, need to clean it up
        String argName = "--founderID";
        String[] goodArguments = new String[]{"s1", "s2",  "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10"};

        CommandLineParser clp1 = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(Collections.singletonList(new InbreedingCoeff(Collections.emptySet())), null)),
                Collections.emptySet());
        CommandLineParser clp2 = new CommandLineArgumentParser(
                new Object(),
                // Adding a default value which should result in different annotations
                Collections.singletonList(new GATKAnnotationPluginDescriptor(Collections.singletonList(new InbreedingCoeff(new HashSet<>(Arrays.asList(goodArguments)))), null)),
                Collections.emptySet());

        List<String> args1 = Stream.of("--annotation", InbreedingCoeff.class.getSimpleName()).collect(Collectors.toList());
        List<String> args2 = Stream.of("--annotation", InbreedingCoeff.class.getSimpleName(), argName, "s1").collect(Collectors.toList());
        Arrays.asList(goodArguments).forEach(arg -> {args1.addAll(Arrays.asList(argName, arg));});

        clp1.parseArguments(nullMessageStream, args1.toArray(new String[args1.size()]));
        List<Annotation> goodOverridingBad = instantiateFilter(clp1);
        clp2.parseArguments(nullMessageStream, args2.toArray(new String[args2.size()]));
        List<Annotation> badOverridingGood = instantiateFilter(clp2);

        Assert.assertEquals(goodOverridingBad.size(), 1);
        Assert.assertEquals(badOverridingGood.size(), 1);
        // assert that the overriding worked in both cases
        Assert.assertEquals(Double.valueOf((String) ((InbreedingCoeff) goodOverridingBad.get(0))
                        .annotate(null, inbreedingCoefficientVC, null)
                        .get(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY)),
                -0.3333333, 0.001, "InbreedingCoefficientScores");
        Assert.assertEquals(((InbreedingCoeff) badOverridingGood.get(0)).annotate(null, inbreedingCoefficientVC, null), Collections.emptyMap());
    }

    //TODO This test is intended to show that the disabling and reenabling an annotation will restore global defaults but unfortunately
    //     this is not the case due to an open issue having to do with instantiating only one annotation object.
    //     See https://github.com/broadinstitute/gatk/issues/3848 for the same issue in read filters.
    @Test (enabled = false)
    public void testDisableDefaultsAndReplaceOnCommandLine() {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(Arrays.asList(new InbreedingCoeff(Sets.newSet("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10")),
                        new DepthPerSampleHC()), Collections.singletonList(StandardAnnotation.class.getSimpleName()))),
                Collections.emptySet());
        List<String> args = Stream.of("-G", StandardAnnotation.class.getSimpleName(), "--"+ StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_ANNOTATIONS).collect(Collectors.toList());

        clp.parseArguments(nullMessageStream, args.toArray(new String[args.size()]));
        VariantAnnotatorEngine vae = new VariantAnnotatorEngine(instantiateFilter(clp), null, Collections.emptyList(), false);

        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());
        Assert.assertTrue(vae.getInfoAnnotations().stream().noneMatch(a -> a.getClass().getSimpleName().equals(DepthPerSampleHC.class.getSimpleName())));
        Assert.assertTrue(vae.getInfoAnnotations().stream().anyMatch(a -> a.getClass().getSimpleName().equals(BaseQualityRankSumTest.class.getSimpleName())));

        // Asserting that the default set InbreedingCoeff has been overridden
        Assert.assertEquals(Double.valueOf((String) ( vae.annotateContext(inbreedingCoefficientVC, new FeatureContext(),null, null, s -> true))
                        .getAttribute(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY)),
                -0, 0.001, "InbreedingCoefficientScores");
    }

    @Test
    public void testEmpty(){
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(null, null)),
                Collections.emptySet());
        String[] args = {};
        clp.parseArguments(nullMessageStream, args);
        List<Annotation> annots = instantiateFilter(clp);
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annots, null, Collections.emptyList(), false);
        Assert.assertTrue(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertTrue(vae.getInfoAnnotations().isEmpty());
    }

    @Test
    public void testExclude(){
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(Collections.singletonList(new Coverage()), null)),
                Collections.emptySet());
        String[] args = {"-AX", Coverage.class.getSimpleName()};
        clp.parseArguments(nullMessageStream, args);
        List<Annotation> annots = instantiateFilter(clp);
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annots, null, Collections.emptyList(), false);
        Assert.assertTrue(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertTrue(vae.getInfoAnnotations().isEmpty());
    }

    @Test
    public void testIncludeDefaultGroupExcludeIndividual(){
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(null, Collections.singletonList(StandardAnnotation.class.getSimpleName()))),
                Collections.emptySet());
        String[] args = {"--annotationsToExclude", Coverage.class.getSimpleName()};
        clp.parseArguments(nullMessageStream, args);
        List<Annotation> annots = instantiateFilter(clp);
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annots, null, Collections.emptyList(), false);

        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());
        //check that Coverage is out
        Assert.assertTrue(vae.getInfoAnnotations().stream().noneMatch(a -> a.getClass().getSimpleName().equals(Coverage.class.getSimpleName())));
    }

    @Test
    public void testIncludeUserGroupExcludeIndividual(){
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(null, null)),
                Collections.emptySet());
        String[] args = {"--annotationsToExclude", Coverage.class.getSimpleName(), "-G", StandardAnnotation.class.getSimpleName()};
        clp.parseArguments(nullMessageStream, args);
        List<Annotation> annots = instantiateFilter(clp);
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annots, null, Collections.emptyList(), false);

        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());
        //check that Coverage is out
        Assert.assertTrue(vae.getInfoAnnotations().stream().noneMatch(a -> a.getClass().getSimpleName().equals(Coverage.class.getSimpleName())));
    }

    //TODO if I decide you need to override then flip this test to assert the opposite
    @Test
    public void testAnnotationGroupOverriding(){
        //TODO this will fail currently and its gross to deal with how group overriding should work
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(null, Collections.singletonList(StandardHCAnnotation.class.getSimpleName()))),
                Collections.emptySet());
        String[] args = {"-G", StandardAnnotation.class.getSimpleName()};
        clp.parseArguments(nullMessageStream, args);
        List<Annotation> annots = instantiateFilter(clp);
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annots, null, Collections.emptyList(), false);

        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());

        //check that Coverage is in and ClippingRankSum is out
        Assert.assertTrue(vae.getInfoAnnotations().stream().anyMatch(a -> a.getClass().getSimpleName().equals(Coverage.class.getSimpleName())));
        Assert.assertTrue(vae.getInfoAnnotations().stream().anyMatch(a -> a.getClass().getSimpleName().equals(ClippingRankSumTest.class.getSimpleName())));
    }
}