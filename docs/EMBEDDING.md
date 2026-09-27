# Embed the PCR engine

Sequence Workbench (Gene-research) is the upper research workflow application.
PrimerBLAST OSS provides local PCR design and specificity, independent of GUI,
HTTP, model provider or user account. Workbench imports the Python library directly.

Use `primerblast_oss.workflows.execute("design", params)` for the same evidence
returned by the agent CLI. `execute("check", params)` screens existing primers.
Database creation is denied by default. Preserve search completeness and coordinate
conventions in client displays; `ok` is not scientific or Wet validation.

`primerblast_oss.workbench_design` preserves the existing Workbench Primer3 settings,
1-based candidate output and exception names. It delegates execution to the same
`design.run_boulder` runner as the native engine. Its code originated in the MIT
Gene-research repository; the original license is included alongside this package.
Different application presets remain explicit, not silently unified.

The standalone GUI defaults to PCR design/check and database preparation.
Legacy sequencing, CAPS, marker and enzyme views remain available through the
legacy-tools switch for compatibility. New broad research workflow UI belongs in
Workbench, not a second independent workbench here. CLI scientific capabilities
remain available to embedded callers.
