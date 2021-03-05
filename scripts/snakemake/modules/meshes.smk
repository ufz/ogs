rule generate_square_mesh:
    output:
        "{mesh_name_prefix,\w+}_{size,\d+}_{lx,\d+}x{ly,\d+}_{type}.vtu"
    shell:
        """
        generateStructuredMesh -e {wildcards.type} \
            --lx {wildcards.lx} --ly {wildcards.ly} \
            --nx {wildcards.size} --ny {wildcards.size} \
            -o {output}
        """
