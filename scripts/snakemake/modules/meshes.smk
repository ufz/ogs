rule generate_square_mesh:
    output:
        "square_{size}_{lx}x{ly}_{type}.vtu"
    shell:
        """
        generateStructuredMesh -e {wildcards.type} \
            --lx {wildcards.lx} --ly {wildcards.ly} \
            --nx {wildcards.size} --ny {wildcards.size} \
            -o {output}
        """
