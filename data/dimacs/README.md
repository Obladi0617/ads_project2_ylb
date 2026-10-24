DIMACS Challenge 9 USA distance graphs

Local inputs are standard DIMACS .gr files for NY, COL, and FLA. They contain directed arcs with integer distance weights.

| Graph | Nodes | Arcs |
|---|---:|---:|
| NY | 264,346 | 733,846 |
| COL | 435,666 | 1,057,066 |
| FLA | 1,070,376 | 2,712,798 |

Source dataset: https://www.diag.uniroma1.it/challenge9/download.shtml
Original download pattern: https://www.diag.uniroma1.it/challenge9/data/USA-road-d/USA-road-d.<NAME>.gr.gz

The original host was unreachable from this environment. The local .gr files were reconstructed from the FROM, TO, and W arrays in https://github.com/sengorajkumar/gpu_graph_algorithms (parser/parser.go reads DIMACS arcs and writes these parallel arrays). The arrays were checked against the published node and arc counts.

The reconstructed .gr files are included in this repository for development convenience. For the course report/submission package, follow the assignment note and list the dataset download links instead of attaching the raw data.