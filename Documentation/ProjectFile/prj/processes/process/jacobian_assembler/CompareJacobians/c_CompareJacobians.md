A Jacobian assembler that assembles the Jacobian in two different ways, compares
the resulting local Jacobians and writes extensive logs in the form of a Python
script if the provided tolerances are exceeded.

Logging (and optionally program termination) is triggered only if both the
absolute and the relative tolerance are exceeded.

# Configuration example

Code snippet:

\code{.xml}
<OpenGeoSysProject>
    <processes>
        <process>
            <jacobian_assembler>
                <type>CompareJacobians</type>
                <jacobian_assembler>
                    <type>Analytical</type>
                </jacobian_assembler>
                <reference_jacobian_assembler>
                    <type>CentralDifferences</type>
                    <component_magnitudes>1 1</component_magnitudes>
                    <relative_epsilons>1e-6 1e-6</relative_epsilons>
                </reference_jacobian_assembler>
                <abs_tol>1e-18</abs_tol>
                <rel_tol>1e-7</rel_tol>
                <fail_on_error>true</fail_on_error>
                <log_file>/tmp/test.log</log_file>
            </jacobian_assembler>
            <!-- ... -->
        </process>
    </processes>
</OpenGeoSysProject>
\endcode
