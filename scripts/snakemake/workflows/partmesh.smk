from snakemake.utils import min_version
min_version("7.3")

module partmesh:
    snakefile: gitlab("bilke/snakemake-partmesh", path="workflow/Snakefile", tag="bbcdc9a721fcda7a1e9619895dce79dfc2b2c905", host="gitlab.opengeosys.org")
    config: config

use rule * from partmesh
