Opens a sequence of mesh tags. Each mesh tag specifies a mesh that have to be
contained in the global mesh list. Output will be generated for each of the
specified meshes. If this list is empty or not present, only the bulk mesh will
be written.

For each mesh a PVD file is created using the naming template specified in the
prefix tag. In this file the VTU files are referenced. These are created
according to the scheme prefix tag suffix tag.
