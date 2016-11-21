The OGS6 input file parameters are documented in the page hierarchy rooted here.

Depending on the type of the parameters the corresponding page titles have
different prefixes, namely:

 - <b>[tag] </b> if the parameter is an XML tag in the input file<br>
   Example: \ref ogs_file_param__prj__linear_solvers
 - <b>[attr]</b> if the parameter is an attribute of an XML tag in the input
   file<br>
   Example: \ref ogs_file_attr__gml__points__point__x
 - <b>[case]</b> either on the top level of the documentation tree, or to
   distinguish between different cases (i.e., sets of configuration options,
   which usually also means different underlying C++ types) in the input
   file.<br>
   Usually one can choose one of several cases by specifying the associated
   <tt>&lt;type&gt;</tt> tag in the input file.<br>
   Example: \ref ogs_file_param__prj__linear_solvers__linear_solver
   and \ref ogs_file_param__parameter__Constant (<tt>&lt;type&gt;Constant&lt;/type&gt;</tt>)
   vs. \ref ogs_file_param__parameter__MeshNode (<tt>&lt;type&gt;MeshNode&lt;/type&gt;</tt>)

The input file parameters are documented within a tree structure (cf. the
navigation tree on the left). This structure resembles the XML document tree of
the input files, however, there are some small differences (see below).

Currently two different input files are documented:
The project file (\ref ogs_file_param__prj "prj")
and the geometry file (\ref ogs_file_param__gml "gml").
All other cases on the top level do not represent separate input files but
rather are shortcuts to certain project file sections, which are provided only in
order to keep the documentation tree flat, they have no other special meaning.

A path in the documentation tree corresponds to a path in the prj or gml input
file, e.g., \ref ogs_file_param__prj__process_variables__process_variable__boundary_conditions__boundary_condition
corresponds to the path <tt>process_variables.process_variable.boundary_conditions.boundary_condition</tt>
in the project file (here parent and child tags are separated by a dot).<br>
Note: The top level XML tags (i.e., <tt>&lt;OpenGeoSysProject&gt;</tt> and <tt>&lt;OpenGeoSysGLI&gt;</tt>)
have been omitted in the entire documentation for the sake of brevity.

There are two exceptions to this rule, both related to the <em>[case]</em>
prefix:
 1. Cases at the top level do not translate to XML tags. The cases
    \ref ogs_file_param__prj and \ref ogs_file_param__gml mark the root of the
    prj and gml input file, respectively.
    All other top level cases, however, belong somewhere in the project file XML tree.
    You can see where they belong in the <em>Additional Info</em> section of the [tag] and
    [attr] pages, e.g., \ref ogs_file_param__boundary_condition__geometry of the
    boundary condition has the full XML tag path
    <tt>process_variables.process_variable.boundary_conditions.boundary_condition.geometry</tt>.
 2. Cases that distinguish types do not translate to XML tags either.
    They enclose a set of configuration settings that can be used at that specific point.<br>
    Example: There is a \ref ogs_file_param__parameter__MeshNode "MeshNode"
    initial condition, which has a \ref
    ogs_file_param__parameter__MeshNode__field_name "field_name" tag.  I.e., if you
    configure the MeshNode parameter, you can (or must) specify the XML tag with the
    path <tt>parameters.parameter.field_name</tt> in the project file.<br> A
    MeshNode parameter can be defined by setting the \ref
    ogs_file_param__parameter__type "parameters.parameter.type" tag to
    <tt>MeshNode</tt>.

# Further Information

 - \subpage ogs_file_param
