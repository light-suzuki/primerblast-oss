# CLI-first application boundary

```text
Local agent / shell -> CLI JSON adapter --+
                                         +-> workflows.py -> existing scientific engines
Browser GUI -> HTTP / job adapter --------+
```

`workflows.py` owns experiment orchestration. `webapp/server.py` owns only HTTP,
job lifecycle, static assets and environment discovery. Compatibility imports for
earlier integrations remain in the server module. New engine orchestration belongs
in the application layer, never in browser code or HTTP handlers. No model provider
or account is required: an agent supplies structured input to the local CLI.

Existing native commands (`design`, `check`, `sequence`, `assay`, etc.) continue
to work. The canonical GUI-compatible JSON entry points are:

```sh
python -m primerblast_oss agent schema
python -m primerblast_oss agent run --input request.json
```

Stdin is supported by `--input -` (the default). Example using only synthetic
primers and a user-supplied local database:

```json
{"request_id":"example-1","operation":"check","params":{"db":["local-reference"],"forward":"ACGTACGTACGTACGTACGT","reverse":"TGCATGCATGCATGCATGCA"}}
```

Supported operations are `design`, `check`, `tile`, `sequence`, `assay`, `markers`,
and `makedb`, sharing the GUI parameter objects and execution functions. `schema`
lists operations and primary input fields; additional existing scientific knobs
remain supported and are documented by the native command's `--help`.
Numeric JSON inputs should be numbers, booleans should be booleans, and `db` and
`primers` should be arrays. Local reference paths and explicit DB/FASTA associations
remain local inputs, never permission to upload their contents.

Stdout contains one UTF-8 JSON envelope with `schema_version`, `ok`, and `result`
or `error`; `request_id` is echoed when supplied. Diagnostics use stderr. Exit code
0 means execution completed, 2 means invalid input, and 1 means execution failure.
An incomplete search remains incomplete inside the result; `ok` never means that
the experiment is validated or that all off-targets were excluded. Existing
coordinates, reports and evidence fields are preserved.

Database creation is blocked by default for agent requests. Only the explicit
`--allow-db-write` option enables `makedb`. No database creation is performed by
discovery, schema inspection or ordinary analysis. The GUI's existing explicit
database-creation action remains available. No new GUI tabs are introduced.
