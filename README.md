# ADS Project 2

## 目录

- src/1.c：Dijkstra 与二叉、三叉、Fibonacci、Binomial 四种优先队列实现。
- data/dimacs/：NY、COL、FLA 道路网络输入文件及来源说明。
- scripts/plot_benchmark.py：从计时 CSV 绘制运行时间图。
- results/：计时、查询、操作计数 CSV，以及生成的 PNG 图。
- build/：编译生成的可执行文件。
- docs/：实验记录和基准分析。

## 编译、运行和画图

请从项目根目录执行：

    g++ -std=c++17 -O2 .\src\1.c -o .\build\dijkstra_benchmark.exe
    .\build\dijkstra_benchmark.exe
    python .\scripts\plot_benchmark.py

.\build\dijkstra_benchmark.exe --count-heap-ops 会单独运行堆内部计数，并更新 results/heap_internal_stats.csv。
The raw .gr graph files are not committed, as required by the assignment. See data/dimacs/README.md for the official download links and expected filenames.