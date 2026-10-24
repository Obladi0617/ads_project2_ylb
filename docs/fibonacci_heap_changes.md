# Fibonacci 堆改动说明

- 日期：2026-10-24
- 对应实现：`src/1.c`
- 改动范围：`fib_consolidate`、`fib_heap_init`，以及 Fibonacci 路径上的计数宏和函数派发

---

## 1. 改动概览

| # | 改动 | 位置 | 目的 |
|---|---|---|---|
| 1 | `consolidate` 不再整表 `memset`，只扫描实际用到的 `[0, hi]` | `fib_consolidate` | 去掉每次 extract-min 的 O(D) 空转 |
| 2 | `degree_table` 改用 `calloc` 初始化 | `fib_heap_init` | 支撑改动 1 的「表在两次整理之间始终为空」不变量 |
| 3 | 内部计数宏改为编译期开关 | `METRIC_ADD` / `METRIC_MAX` | 计时版不再执行 27 亿次计数分支 |
| 4 | Dijkstra 对 Fibonacci 堆宏特化 | `DEFINE_DIJKSTRA_QUERY` | 堆操作内联，去掉每次操作的函数指针派发 |

其中 1、2 是 Fibonacci 专属改动；3、4 影响四种堆，这里只记录它们在 Fibonacci 路径上的效果。

---

## 2. 原实现的问题

### 2.1 每次 extract-min 都清空并扫描整张 degree 表

原代码：

```c
memset(heap->degree_table,0,
       (size_t)heap->degree_capacity*sizeof(FibNode*));

/* ... 链接同 degree 的树 ... */

for(int i=0;i<heap->degree_capacity;i++){
    FibNode *x=heap->degree_table[i];
    if(x==NULL) continue;
    /* 把 x 放回根表 */
}
```

`degree_capacity` 是按**整张图的点数**算的：

```c
int capacity=3;
for(int x=n;x>1;x>>=1){
    capacity+=2;
}
```

| 图 | 点数 | `degree_capacity` | extract-min 总次数 | 全表扫描迭代 | memset 字节量 |
|---|---:|---:|---:|---:|---:|
| NY | 264,346 | 39 | 135,502,866 | 5.28 × 10⁹ | 42.3 GB |
| COL | 435,666 | 39 | 220,809,144 | 8.61 × 10⁹ | 68.9 GB |
| FLA | 1,070,376 | 43 | 531,602,948 | 22.86 × 10⁹ | 182.9 GB |

但实际每次整理只有约 **12 个根**（NY：`root_list_steps` 1.686 × 10⁹ / 135.5 × 10⁶ ≈ 12.4），
最大 degree 约 9–10。也就是说，每次 extract-min 有 **2/3 以上的表扫描是纯空转**。

### 2.2 表项在两次整理之间没有「保持为空」这个约定

原实现靠每次开头的 `memset` 保证表是空的，所以不能在中间省掉它。
要把扫描范围缩到 `[0, hi]`，必须先建立一条更强的不变量（见 3.3）。

### 2.3 计数分支在 Fibonacci 热路径里

`fib_node_less`、`fib_add_root`、`fib_link`、`fib_heap_extract_min` 里的
`METRIC_ADD` 原本是**运行时**判断：

```c
#define METRIC_ADD(field,value) do{ \
    if(active_heap_metrics!=NULL){ \
        active_heap_metrics->field+=(value); \
    } \
}while(0)
```

NY 上 Fibonacci 的比较次数是 **2.745 × 10⁹**，这些全局指针检查和分支都在计时运行里执行。

### 2.4 每次堆操作都走函数指针

原 `dijkstra_query` 通过 `MinPQ` 的函数指针表调用
`empty / insert / extract_min / decrease_key`。NY 上 Fibonacci 路径大约有
2.7 亿次 `empty` 检查、1.36 亿次 insert、1.36 亿次 extract、0.19 亿次 decrease-key，
全部无法内联。

---

## 3. 改动一：consolidate 只扫描 `[0, hi]`

### 3.1 改后代码

```c
/* The degree table stays all-NULL between consolidations: entries
   consumed during linking are cleared immediately, and surviving
   entries are cleared below as they move back to the root list.
   Only the touched range [0,hi] is scanned, which removes the former
   O(degree_capacity) memset and full-table scan per extract-min. */
int hi=-1;

for(int i=0;i<root_count;i++){
    FibNode *x=heap->root_buffer[i];
    if(x->parent!=NULL){
        continue;
    }

    int degree=x->degree;
    while(degree<heap->degree_capacity&&
          heap->degree_table[degree]!=NULL){
        FibNode *y=heap->degree_table[degree];

        if(fib_node_less(y,x)){
            FibNode *temp=x;
            x=y;
            y=temp;
        }

        heap->degree_table[degree]=NULL;   /* 被消费的项立即清掉 */
        fib_link(y,x);
        degree=x->degree;
    }

    if(degree>=heap->degree_capacity){
        fprintf(stderr,"Fibonacci heap degree bound exceeded\n");
        return 0;
    }

    heap->degree_table[degree]=x;
    if(degree>hi){
        hi=degree;
    }
}

heap->min=NULL;

for(int i=0;i<=hi;i++){
    FibNode *x=heap->degree_table[i];
    if(x==NULL){
        continue;
    }

    heap->degree_table[i]=NULL;            /* 存活项也清掉 */
    x->left=x;
    x->right=x;
    x->parent=NULL;
    x->mark=0;
    fib_add_root(heap,x);
}
```

### 3.2 改动点

- 删除每次调用开头的 `memset`。
- 链接循环里，被消费的表项**当场置 NULL**。
- 用一个 `hi` 记录本次写入过的最大下标，最后只扫描 `0..hi`。
- 扫描时把存活项也**当场置 NULL**。

### 3.3 正确性：表在两次整理之间始终为空

需要证明：`fib_consolidate` 正常返回时，`degree_table` 的所有项都是 `NULL`。

1. **初始状态**：`fib_heap_init` 用 `calloc` 分配 `degree_table`，全为 `NULL`。
2. **归纳假设**：上一次整理返回后，表全为 `NULL`。
3. **本次整理**：
   - 链接循环中，每个被读到并消费的 `table[degree]` 都执行 `table[degree]=NULL`；
   - 每个新写入的 `table[degree]=x` 都有 `degree ≤ hi`；
   - 收尾循环遍历 `0..hi`，对每个非空项执行 `table[i]=NULL`，并把它放回根表。
4. 因此本次所有写过的下标都已被清空；结合归纳假设，返回时表全为 `NULL`。

两个边界情况：

- `hi = -1`（本次没有任何根被写入）：收尾循环不执行，表保持全空。
  正常堆里所有根都有 `parent == NULL`，这种情况只可能出现在非法输入上。
- 度数越界走 `return 0`：此时表里可能残留项。但该错误会让整个 benchmark
  立即终止，堆随后被销毁，不会再次整理，所以不影响正确性。

### 3.4 复杂度变化

| | 改前（每次 extract-min） | 改后 |
|---|---|---|
| 清零 | `O(degree_capacity)` = 39 / 39 / 43 次写 | 0 |
| 收尾扫描 | `O(degree_capacity)` | `O(hi+1)`，实测约 11 |
| 链接 | `O(root_count)` | `O(root_count)`（不变） |

按 NY 的内层数据估算，收尾扫描从 39 项降到约 11 项，
FLA 从 43 项降到约 11 项，**单次 extract-min 的固定开销减少约 70%**。

---

## 4. 改动二：`degree_table` 用 `calloc`

```c
/* 改前 */
heap->degree_table=(FibNode**)malloc((size_t)capacity*sizeof(FibNode*));

/* 改后 */
heap->degree_table=(FibNode**)calloc((size_t)capacity,sizeof(FibNode*));
```

这是改动一的前提：没有开头的 `memset` 后，必须保证初次分配时表是全空的。

---

## 5. 改动三：计数宏改为编译期开关

```c
#ifdef COUNT_HEAP_OPS
static HeapInternalStats *active_heap_metrics=NULL;
#define METRIC_ADD(field,value) do{ \
    if(active_heap_metrics!=NULL){ \
        active_heap_metrics->field+=(value); \
    } \
}while(0)
#define METRIC_MAX(field,value) do{ \
    if(active_heap_metrics!=NULL&& \
       active_heap_metrics->field<(long long)(value)){ \
        active_heap_metrics->field=(long long)(value); \
    } \
}while(0)
#else
#define METRIC_ADD(field,value) do{}while(0)
#define METRIC_MAX(field,value) do{}while(0)
#endif
```

计时版用 `-O2` 编译，`METRIC_ADD` / `METRIC_MAX` 直接展开为空，
`fib_node_less`、`fib_add_root`、`fib_link`、`extract_min` 里不再有指针检查和加法。

计数版用 `-O2 -DCOUNT_HEAP_OPS` 编译，`--count-heap-ops` 也只在这个版本里可用；
用计时版执行该参数会提示重新编译。

---

## 6. 改动四：Fibonacci 的堆操作改为宏特化

`dijkstra_query` 现在由宏生成四份，Fibonacci 一份的实例化是：

```c
DEFINE_DIJKSTRA_QUERY(fibonacci,FibHeap,
                      fib_heap_empty,fib_heap_insert,
                      fib_heap_extract_min,fib_heap_decrease_key,
                      fib_pq_clear)
```

展开后，`fib_heap_empty / fib_heap_insert / fib_heap_extract_min /
fib_heap_decrease_key` 都是编译期确定的直接调用，编译器可以内联；
不再经过 `MinPQ` 的函数指针表。原 `fib_pq_decrease_key` 包装函数随之删除。

---

## 7. 这次没有改、但可以继续做的

这几项没有动，保持本次改动可验证、可回退：

1. **节点体积**：`FibNode` 目前 72 字节（`long long key` + 3 个 int + 4 个指针）。
   改成 32 位 key、位域 `mark/in_heap`、32 位下标代替指针后，可以压到 32–40 字节，
   对 1M 节点的指针追逐会有明显收益。
2. **距离类型**：`dist[]`、`Item.key`、`FibNode.key` 都是 64 位。
   USA 道路网络的最长最短路远小于 2³¹，改成 32 位可以减半内存带宽。
3. **延迟整理（lazy consolidation）**：当根表长度小于 `2*log2(size)` 时可以跳过
   consolidate，`min` 由 `fib_add_root` 维护；在 Dijkstra 这种小队列场景可能进一步减少整理次数。
4. **degree 表容量**：仍按整图 `n` 计算（39/43），可以改成按当前堆大小动态计算，
   进一步降低 `hi` 的上界；本次改动的收益已经不依赖这一点。
5. **`root_buffer`**：仍按 `n` 个指针分配（FLA 约 8.6 MB），可以改成按当前堆大小增长。

---

## 8. 验证

| 测试 | 结果 |
|---|---|
| 50 查询冒烟（NY，四种堆） | 50/50 答案一致 |
| 1000 查询（NY，四种堆） | 1000/1000 答案一致 |
| 三图全量 | 见下一节 |

冒烟命令（把查询数临时改成 50、只跑 NY）：

```powershell
# 编译计时版
g++ -std=c++17 -O2 .\src\1.c -o .\build\dijkstra_benchmark.exe

# 编译计数版
g++ -std=c++17 -O2 -DCOUNT_HEAP_OPS .\src\1.c -o .\build\dijkstra_benchmark_count.exe

# 计时
.\build\dijkstra_benchmark.exe

# 只统计堆内部操作
.\build\dijkstra_benchmark_count.exe --count-heap-ops
```

---

## 9. 实测结果

### 9.1 本机基线（改动前，NY，1000 查询）

| 堆 | 时间 |
|---|---:|
| Binary | 27.020 s |
| Ternary | 28.443 s |
| Fibonacci | 68.906 s |
| Binomial | 43.290 s |

环境：Windows 11，MinGW g++ 8.1.0，`-O2`，单线程。
同一二进制连续两轮的漂移约 4–8%，因此小于该量级的差异不应从单轮结果下结论。

### 9.2 优化版（NY / COL / FLA，1000 查询）

> 全量运行完成后补充三张图的实际时间、相对基线的变化，以及对应的
> `results/benchmark_results.csv` 和 `results/benchmark_runtime.png`。

---

## 10. 复现与回退

- 原始版本保存在 `ads_project2_ylb_scratch/orig.c`（工作副本，未提交）。
- 改动集中在 `src/1.c` 的 `fib_consolidate`、`fib_heap_init`、
  `METRIC_ADD/METRIC_MAX` 和 `DEFINE_DIJKSTRA_QUERY`。
- 回退时只需恢复这四个位置；其余堆的实现未因 Fibonacci 改动而改变。
