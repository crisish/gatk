package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import com.google.common.annotations.VisibleForTesting;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.utils.Utils;
import org.reflections.ReflectionUtils;

import java.lang.reflect.Modifier;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * A plugin descriptor for managing the dynamic discovery of both @{InfoFieldAnnotation} and @{GenotypeAnnotation} objects
 * within the @{org.broadinstitute.hellbender.tools.walkers.annotator} package which handles integrating annotation specific
 * arguments from the command line with tool specified defaults.
 *
 * Unlike @{GATKReadFilterPluginDescriptor} annotation order is not important and thus argument order is not guarinteed to be
 * preserved in all cases, especially when group annotations are involved.
 */
public class GATKAnnotationPluginDescriptor  extends CommandLinePluginDescriptor<Annotation> {
    private static final String pluginPackageName = "org.broadinstitute.hellbender.tools.walkers.annotator";
    private static final Class<?> pluginBaseClass = org.broadinstitute.hellbender.tools.walkers.annotator.Annotation.class;

    protected transient Logger logger = LogManager.getLogger(this.getClass());

    @VisibleForTesting
    @ArgumentCollection
    private final VariantAnnotationArgumentCollection userArgs;

    // Map of Annotation (simple) class names to the corresponding discovered plugin instance
    private final Map<String, Annotation> allDiscoveredAnnotations = new HashMap<>();

    // Map of Annotation (simple) class names to the corresponding default plugin instance
    // it is a LinkedHashMap because we want to remember the order in which these were provided, and also keep the
    // actual instances in case they have any additional state provided by the tool
    // when they were created
    private final Map<String, Annotation> toolDefaultAnnotations = new LinkedHashMap<>();
    private final Set<String> toolDefaultGroups = new HashSet<>();

    // Set of dependent args for which we've seen values (requires predecessor)
    private final Set<String> requiredPredecessors = new HashSet<>();

    private final Map<String, List<Annotation>> discoveredGroups = new HashMap<>();

    /**
     * @param userArgs           Argument collection to control the exposure of the command line arguments.
     * @param toolDefaultFilters Default filters that may be supplied with arguments
     *                           on the command line. May be null.
     */
    public GATKAnnotationPluginDescriptor(final VariantAnnotationArgumentCollection userArgs, final List<Annotation> toolDefaultFilters, final List<String> toolDefaultGroups) {
        this.userArgs = userArgs;
        if (null != toolDefaultFilters) {
            toolDefaultFilters.forEach(f -> {
                final Class<? extends Annotation> rfClass = f.getClass();
                // anonymous classes have a 0-length simple name, and thus cannot be accessed or
                // controlled by the user via the command line, but they should still be valid
                // as default filters, so use the full name to ensure that their map entries
                // don't clobber each other
                String className = rfClass.getSimpleName();
                if (className.length() == 0) {//TODO possibly unecissary for annotations
                    className = rfClass.getName();
                }
                toolDefaultAnnotations.put(className, f);
            });
        }
        if (null != toolDefaultGroups) {
            this.toolDefaultGroups.addAll(toolDefaultGroups);
        }
    }

    /**
     * @param toolDefaultFilters Default filters that may be supplied with arguments
     *                           on the command line. May be null.
     */
    public GATKAnnotationPluginDescriptor(final List<Annotation> toolDefaultFilters, final List<String> toolDefaultGroups) {
        this(new VariantAnnotationArgumentCollection(Collections.emptyList(), Collections.emptyList(), Collections.emptyList()), toolDefaultFilters, toolDefaultGroups);
    }

    @Override
    public Predicate<Class<?>> getClassFilter() {
        return c -> {
            // don't use the Annotation base class, it's inner classes, or the unit tests
            return !c.getName().equals(this.getPluginClass().getName()) &&
                    !Modifier.isAbstract(c.getModifiers()) &&
                    !c.getName().contains("UnitTest$");
        };
    }

    /**
     * Return a display name to identify this plugin to the user
     *
     * @return A short user-friendly name for this plugin.
     */
    @Override
    public String getDisplayName() {
        return StandardArgumentDefinitions.ANNOTATION_LONG_NAME;
    }

    /**
     * @return the class object for the base class of all plugins managed by this descriptor
     */
    @Override
    public Class<?> getPluginClass() {
        return pluginBaseClass;
    }

    /**
     * A list of package names which will be searched for plugins managed by the descriptor.
     *
     * @return
     */
    @Override
    public List<String> getPackageNames() {
        return Collections.singletonList(pluginPackageName);
    }

    /**
     * Return an instance of the specified pluggable class. The descriptor should
     * instantiate or otherwise obtain (possibly by having been provided an instance
     * through the descriptor's constructor) an instance of this plugin class.
     * The descriptor should maintain a list of these instances so they can later
     * be retrieved by {@link #getAllInstances}.
     * <p>
     * In addition, implementations should recognize and reject any attempt to instantiate
     * a second instance of a plugin that has the same simple class name as another plugin
     * controlled by this descriptor (which can happen if they have different qualified names
     * within the base package used by the descriptor) since the user has no way to disambiguate
     * these on the command line).
     *
     * @param pluggableClass a plugin class discovered by the command line parser that
     *                       was not rejected by {@link #getClassFilter}
     * @return the instantiated object that will be used by the command line parser
     * as an argument source
     * @throws IllegalAccessException
     * @throws InstantiationException
     */
    //TODO
    @Override
    public Object getInstance(Class<?> pluggableClass) throws IllegalAccessException, InstantiationException {
        Annotation annot = null;
        final String simpleName = pluggableClass.getSimpleName();

        if (allDiscoveredAnnotations.containsKey(simpleName)) {
            // we found a plugin class with a name that collides with an existing class;
            // plugin names must be unique even across packages
            throw new IllegalArgumentException(
                    String.format("A plugin class name collision was detected (%s/%s). " +
                                    "Simple names of plugin classes must be unique across packages.",
                            pluggableClass.getName(),
                            allDiscoveredAnnotations.get(simpleName).getClass().getName())
            );
        } else if (toolDefaultAnnotations.containsKey(simpleName)) {
            // an instance of this class was provided by the tool as one of it's default filters;
            // use the default instance as the target for command line argument values
            // rather than creating a new one, in case it has state provided by the tool
            annot = toolDefaultAnnotations.get(simpleName);
        } else {
            annot = (Annotation) pluggableClass.newInstance();
        }

        // Add all filters to the allDiscoveredReadFilters list, even if the instance came from the
        // tool defaults list (we want the actual instances to be shared to preserve state)
        allDiscoveredAnnotations.put(simpleName, annot);

        Class<?>[] interfaces = annot.getClass().getInterfaces();
        for (Class<?> inter : interfaces) {
            // Following with how groups are currently defined and discovered, namely they are interfaces that
            // extend Annotation. This might not be entirely desirable as it results in ReducibleAnnotation being
            // a definable group even though this is not an entirely helpful list.
            if ((inter != pluginBaseClass) && (pluginBaseClass.isAssignableFrom(inter))) {
                List<Annotation> list = discoveredGroups.containsKey(inter.getSimpleName()) ? discoveredGroups.get(inter.getSimpleName()) : new ArrayList<>();
                list.add(annot);
                discoveredGroups.put(inter.getSimpleName(), list);

                // If its a valid group, check whether the tool requested that group and add it to default annotations
                if (toolDefaultGroups.contains(inter.getSimpleName()) && !toolDefaultAnnotations.containsKey(simpleName)) {
                    toolDefaultAnnotations.put(simpleName, annot);
                }
            }
        }

        return annot;
    }

    /**
     * Return the allowed values for readFilterNames/disableReadFilter for use by the help system.
     *
     * @param longArgName long name of the argument for which help is requested
     * @return
     */
    @Override
    public Set<String> getAllowedValuesForDescriptorArgument(String longArgName) {
        if (longArgName.equals(StandardArgumentDefinitions.ANNOTATION_LONG_NAME)) {
            return allDiscoveredAnnotations.keySet();
        }
        if (longArgName.equals(StandardArgumentDefinitions.ANNOTATIONS_TO_EXCLUDE_LONG_NAME)) {
            return toolDefaultAnnotations.keySet();
        }
        throw new IllegalArgumentException("Allowed values request for unrecognized string argument: " + longArgName);
    }

    @Override
    public boolean isDependentArgumentAllowed(final Class<?> predecessorClass) {
        // Make sure the predecessor for a dependent argument was either specified on the command line or
        // is a tool default, otherwise reject it.
        // NOTE: This method is called by the CLP during parsing at the time the depended argument is seen
        // on the command line. Even if this check passes at the time this method is called, its possible
        // for the user to subsequently disable the required predecessor. That case is caught during final
        // validation done by the validateArguments method.
        String predecessorName = predecessorClass.getSimpleName();
        boolean isAllowed = (userArgs.annotationsToUse.contains(predecessorName))
                || (toolDefaultAnnotations.get(predecessorName) != null);
        if (!isAllowed) {
            // Need to check whether any of the annotations have been included in one of the group dependencies
            // TODO an alternative would be to use ClassUtils.knownSubInterfaceSimpleNames(Annotation.class); to discover the groups ahead of time
            for (String group : Stream.of(userArgs.annotationGroupsToUse, toolDefaultGroups)
                    .flatMap(Collection::stream)
                    .collect(Collectors.toList())) {
                if (discoveredGroups.containsKey(group)
                        && discoveredGroups.get(group).stream().anyMatch(s -> s.getClass().getSimpleName().equals(predecessorName))) {
                    isAllowed = true;
                    break;
                }
            }
        }
        if (isAllowed) {
            // Keep track of the ones we allow so we can validate later that they weren't subsequently disabled
            requiredPredecessors.add(predecessorName);
        }
        return isAllowed;
    }

    /**
     * Validate the list of arguments and reduce the list of read filters to those
     * actually seen on the command line. This is called by the command line parser
     * after all arguments have been parsed. Tries to catch most cases where the user
     * provides potentially confusing input.
     */
    @Override
    public void validateArguments() throws CommandLineException {
        // throw if a tool default annotation group is invalid, this is the tool authors fault and should be a hard failure
        if (!discoveredGroups.keySet().containsAll(toolDefaultGroups)) {
            toolDefaultGroups.removeAll(discoveredGroups.keySet());
            throw new GATKException(String.format("The tool requested annotation group(s) \"%s\" are invalid annotation groups", Utils.join("\", \"", toolDefaultGroups)));
        }

        // throw if a annotation is *enabled* more than once by the user
        final Set<String> duplicateUserEnabledFilterNames = Utils.getDuplicatedItems(userArgs.annotationsToUse);
        if (!duplicateUserEnabledFilterNames.isEmpty()) {
            throw new CommandLineException.BadArgumentValue(
                    String.format("The annotation(s) are enabled more than once: %s",
                            Utils.join(", ", duplicateUserEnabledFilterNames)));
        }

        // throw if a annotation is *disabled* more than once by the user
        final Set<String> duplicateDisabledUserFilterNames = Utils.getDuplicatedItems(userArgs.annotationsToExclude);
        if (!duplicateDisabledUserFilterNames.isEmpty()) {
            throw new CommandLineException.BadArgumentValue(
                    String.format("The annotation(s) are disabled more than once: %s",
                            Utils.join(", ", duplicateDisabledUserFilterNames)));
        }

        // throw if a annotation is both enabled *and* disabled by the user
        final Set<String> enabledAndDisabled = new HashSet<>(userArgs.annotationsToUse);
        enabledAndDisabled.retainAll(userArgs.annotationsToExclude);
        if (!enabledAndDisabled.isEmpty()) {
            final String badFiltersList = Utils.join(", ", enabledAndDisabled);
            throw new CommandLineException(
                    String.format("The annotation(s): %s are both enabled and disabled", badFiltersList));
        }

        // throw if a disabled annotation doesn't exist; warn if it wasn't enabled by the tool in the first place
        userArgs.annotationsToExclude.forEach(s -> {
            if (!allDiscoveredAnnotations.containsKey(s)) {
                throw new CommandLineException.BadArgumentValue(String.format("Disabled filter (%s) does not exist", s));
            } else if (!toolDefaultAnnotations.containsKey(s)) {
                logger.warn(String.format("Disabled annotation (%s) is not enabled by this tool", s));
            }
        });

        // warn if an annotation is both default and enabled by the user
        final Set<String> redundantAnnots = new HashSet<>(toolDefaultAnnotations.keySet());
        redundantAnnots.retainAll(userArgs.annotationsToUse);
        redundantAnnots.forEach(
                s -> {
                    logger.warn(String.format("Redundant enabled annotation (%s) is enabled for this tool by default", s));
                });

        // warn if an annotation is both default and enabled by the user
        final Set<String> redundantGroups = new HashSet<>(toolDefaultGroups);
        redundantGroups.retainAll(userArgs.annotationGroupsToUse);
        redundantGroups.forEach(
                s -> {
                    logger.warn(String.format("Redundant enabled annotation group (%s) is enabled for this tool by default", s));
                });

        // Throw if args were specified for a filter that was also disabled, or that was not enabled by the
        // tool by default.
        //
        // Note that this is also checked during command line argument parsing, but needs to be checked again
        // here. Whenever the command line parser sees a dependent argument on the command line, it delegates
        // back to the descriptor's isDependentArgumentAllowed method to allow it to validate that the predecessor
        // for that dependent argument has been supplied, either by a default read filter, or by an explicitly
        // enabled read filter. However, its possible for the user to subsequently try to disable that
        // predecessor, which is what we want to catch here.
        //
        userArgs.annotationsToExclude.forEach(s -> {
            if (requiredPredecessors.contains(s)) {
                String message = String.format("Values were supplied for (%s) that is also disabled", s);
                if (toolDefaultAnnotations.containsKey(s)) {
                    // NOTE: https://github.com/broadinstitute/barclay/issues/23
                    // This is a special case to work around the issue where we can't really tell if the
                    // predecessor was added as a result of a user-provided value, or a default value. The
                    // CLP doesn't distinguish, so we only warn here for now.
                    logger.warn(message);
                } else {
                    throw new CommandLineException(message);
                }
            }
        });

        // throw if an annotation name was specified that has no corresponding instance
        userArgs.annotationsToUse.forEach(s -> {
            Annotation trf = allDiscoveredAnnotations.get(s);
            if (null == trf) {
                if (!toolDefaultAnnotations.containsKey(s)) {
                    throw new CommandLineException("Unrecognized annotation name: " + s);
                }
            }
        });

        // throw if a annotation name was specified that has no corresponding instance
        userArgs.annotationGroupsToUse.forEach(s -> {
            if (!discoveredGroups.containsKey(s)) {
                throw new CommandLineException("Unrecognized annotation group name: " + s);
            }
        });
    }

    /**
     * Get the list of default plugins used for this instance of this descriptor. Used for help/doc generation.
     * <p>
     * NOTE: this method does not account for disabled default filters and just return ALL default instances.
     * The refactored interface in Barclay changes it's contract to allows returning a list with only 'enabled' default
     * instances. We'll change the implementation when we integrate the updated interface.
     */
    @Override
    public List<Object> getDefaultInstances() {
        return new ArrayList<>(toolDefaultAnnotations.values());
    }

    /**
     * Pass back the list of Annotation instances that were actually seen on the command line in the annotations
     * originating from groups first, followed by individual annotaitons in order. Its possible for this to return
     * a filter that was originally included in the list of tool defaults only in the case where the user also
     * specifies it on the command line.
     * <p>
     * NOTE: this method is somewhat misnamed in that it doesn't return ALL instances since it leaves out
     * default filters (Except as noted above). The refactored interface in Barclay renames this method and
     * changes it's contract. We'll change the implementation when we integrate the updated interface.
     */
    @Override
    public List<Annotation> getAllInstances() {
        final LinkedHashSet<Annotation> annotations = new LinkedHashSet<>();
        userArgs.annotationGroupsToUse.forEach(s -> {
            List<Annotation> as = discoveredGroups.get(s);
            annotations.addAll(as);
        });
        userArgs.annotationsToUse.forEach(s -> {
            Annotation rf = allDiscoveredAnnotations.get(s);
            annotations.add(rf);
        });
        return new ArrayList<>(annotations);
    }


    /**
     * Return the class representing the instance of the plugin specified by {@code pluginName}
     *
     * @param pluginName Name of the plugin requested
     * @return Class object for the plugin instance requested
     */
    @Override
    public Class<?> getClassForInstance(final String pluginName) {
        return allDiscoveredAnnotations.get(pluginName).getClass();
    }

    /**
     * Merge the default filters with the users's command line read filter requests, then initialize
     * the resulting filters. Specifically, unless the user disables all tool default annotations it will
     * first add all the tool enabled annotations which were not individually blocked by the user and then
     * adds in annotations defined by the users specified groups, then individual annotations.
     *
     * @return Single merged read filter.
     */
    public Collection<Annotation> getMergedAnnotations() {
        final SortedSet<Annotation> annotations = new TreeSet<>(Comparator.comparing(t -> t.getClass().getSimpleName()));

        if (!userArgs.disableToolDefaultAnnotaitons) {
            annotations.addAll(toolDefaultAnnotations.values().stream()
                    .filter(t -> !userArgs.annotationsToExclude.contains(t.getClass().getSimpleName()))
                    .collect(Collectors.toList()));
        }
        for (final Annotation annot : allDiscoveredAnnotations.values()) {
            if (!userArgs.annotationsToExclude.contains(annot.getClass().getSimpleName())) {
                //if any group matches requested groups, it's in
                //TODO this reflection should be abstracted away to somewhere nicer
                @SuppressWarnings("unchecked") final Set<Class<?>> annotationGroupsForT = ReflectionUtils.getAllSuperTypes(annot.getClass(), sup -> sup.isInterface() && discoveredGroups.keySet().contains(sup.getSimpleName()));
                if (annotationGroupsForT.stream().anyMatch(iface -> userArgs.annotationGroupsToUse.contains(iface.getSimpleName()))) {
                    if (!annotations.contains(annot)) {
                        annotations.add(annot);
                    }
                }
                // If the user explicitly requests an annotation then we want to override its tool-default configuration
                if (userArgs.annotationsToUse.contains(annot.getClass().getSimpleName())) {
                    annotations.add(annot);
                }
            }
        }

        return annotations;
    }
}
