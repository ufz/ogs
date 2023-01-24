This tag specifies that residuum vectors should be assembled and output for a set of
submeshes of the entire simulation domain.

The same settings/child tags as for \ref ogs_file_param__prj__time_loop__output "&lt;output&gt;" can be used, here.
I.e., this tag works mostly like an additional output section next to an
already existing \ref ogs_file_param__prj__time_loop__output "&lt;output&gt;"
or \ref ogs_file_param__prj__time_loop__outputs "&lt;outputs&gt;"
tag.
See also the examples below.

\attention There is one important difference between this tag and the "regular" output
configuration:
The \ref ogs_file_param__prj__time_loop__output__meshes "&lt;meshes&gt;"
must be non-overlapping and must cover the entire simulation domain.
OGS will refuse to run if these requirements are not fulfilled.

\todo As of January 2023 the submesh residuum output has not been implemented
for all processes, yet. And it has not been tested together with other OGS
features, such as staggered coupling, domain decomposition and domain
deactivation.
