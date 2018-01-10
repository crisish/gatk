package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AnnotatedVariantProducer;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryDataBundle;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervalTree;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.stream.Collectors;

public class ImpreciseVariantDetector {

    public static List<VariantContext> detectImpreciseVariantsAndAddReadAnnotations(final SvDiscoveryDataBundle svDiscoveryDataBundle,
                                                                                    List<VariantContext> annotatedVariants) {
        final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks = svDiscoveryDataBundle.evidenceTargetLinks;

        if (evidenceTargetLinks != null) {
            annotatedVariants = processEvidenceTargetLinks(
                    annotatedVariants,
                    evidenceTargetLinks,
                    svDiscoveryDataBundle.metadata,
                    svDiscoveryDataBundle.referenceBroadcast.getValue(),
                    svDiscoveryDataBundle.discoverStageArgs,
                    svDiscoveryDataBundle.toolLogger);
        }
        return annotatedVariants;
    }

    /**
     * Uses the input EvidenceTargetLinks to either annotate the variants called from assembly discovery with split
     * read and read pair evidence, or to call new imprecise variants if the number of pieces of evidence exceeds
     * a given threshold.
     */
    @VisibleForTesting
    public static List<VariantContext> processEvidenceTargetLinks(final List<VariantContext> assemblyDiscoveredVariants,
                                                                  final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                                                  final ReadMetadata metadata,
                                                                  final ReferenceMultiSource reference,
                                                                  final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection parameters,
                                                                  final Logger localLogger) {
        final int originalEvidenceLinkSize = evidenceTargetLinks.size();
        final List<VariantContext> result = assemblyDiscoveredVariants
                .stream()
                .map(variant -> AnnotatedVariantProducer.annotateWithImpreciseEvidenceLinks(
                        variant,
                        evidenceTargetLinks,
                        reference.getReferenceSequenceDictionary(null),
                        metadata, parameters.assemblyImpreciseEvidenceOverlapUncertainty))
                        .collect(Collectors.toList());
        localLogger.info("Used " + (originalEvidenceLinkSize - evidenceTargetLinks.size()) + " evidence target links to annotate assembled breakpoints");

        final List<VariantContext> impreciseVariants =
                Utils.stream(evidenceTargetLinks)
                        .map(p -> p._2)
                        .filter(EvidenceTargetLink::isImpreciseDeletion)
                        .filter(e -> e.getReadPairs() + e.getSplitReads() > parameters.impreciseEvidenceVariantCallingThreshold)
                        .map(e -> createImpreciseDeletionVariant(e, metadata, reference))
                        .collect(Collectors.toList());

        localLogger.info("Called " + impreciseVariants.size() + " imprecise deletion variants");
        result.addAll(impreciseVariants);
        return result;
    }

    private static VariantContext createImpreciseDeletionVariant(final EvidenceTargetLink e,
                                                                 final ReadMetadata metadata,
                                                                 final ReferenceMultiSource reference) {
        final SvType svType = new SimpleSVType.ImpreciseDeletion(e, metadata);
        return AnnotatedVariantProducer
                .produceAnnotatedVcFromEvidenceTargetLink(e, svType, metadata, reference);
    }
}
