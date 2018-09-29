Required tools:

- `time`
- `diff`
- `valgrind`
- `grep`
- `bash`

Test data is searched for in `ogs6-sources/../ogs6-data`.

In the build directory run `ctest` with `–-output-on-failure` and `-R` for includes and `-E` for excludes:

```bash
ctest --output-on-failure -R "MEMCHECK|VALGRIND" -E DIFF
```

Wrapper and tester are implemented in `AddTest.cmake`.
