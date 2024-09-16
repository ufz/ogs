The "id" is an optional parameter.
It starts with 0, and is referring to different BHEs based on their "MaterialID".
Then those BHEs will follow the configurations defined in this block.
A mixed definition with and without id is not supported.
If "id" is not given, then it will be set automatically with an ascending order, also starting from 0.
With the "id" present, multiple BHEs can share the same configurations defined in a single block without repeating it.
To specify the same definition for all BHE in the project, the notation id="*" can be used.
