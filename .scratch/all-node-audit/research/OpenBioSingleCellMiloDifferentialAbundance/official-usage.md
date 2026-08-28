# OpenBioSingleCellMiloDifferentialAbundance：官方用法核对

调研日期：2026-08-28  
核对目标：Pertpy `Milo` 1.3.0、原始 Milo/miloR、Scanpy 1.12 系列与 edgeR 准似然工作流。本文只记录改前证据，不代表当前节点已经满足这些要求。

## 1. 版本与一手资料

本次只采用官方文档、带版本标签的官方源码、包官方页面和方法论文。

| 主题 | 核验来源 | 本次使用的事实 |
| --- | --- | --- |
| Pertpy 当前发行版 | [PyPI: pertpy](https://pypi.org/project/pertpy/)、[Pertpy 1.3.0 release](https://github.com/scverse/pertpy/releases/tag/1.3.0) | 2026-08-24 发布 1.3.0；本记录以 `1.3.0` 标签源码为可复核基准，而不是浮动的 `main`。 |
| Pertpy Milo API | [官方 API](https://pertpy.readthedocs.io/en/latest/api/tools/pertpy.tools.Milo.html)、[1.3.0 `_milo.py` 源码](https://github.com/scverse/pertpy/blob/1.3.0/src/pertpy/tools/_milo.py) | 精确签名、默认求解器、读写槽位、结果列及空间 FDR 实现。 |
| Pertpy Milo 教程 | [官方教程](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/milo.html) | 推荐执行顺序、Sample 计数、表示空间和邻居参数的示例、edgeR/PyDESeq2 差异。 |
| Pertpy 安装与变更 | [安装说明](https://pertpy.readthedocs.io/en/stable/installation.html)、[changelog](https://pertpy.readthedocs.io/en/latest/changelog.html) | edgeR 路径的 R/rpy2 依赖；1.0–1.3 间求解器、混合模型及残余自由度行为发生过变化。 |
| Pertpy 论文 | Heumos et al., [*Pertpy: an end-to-end framework for perturbation analysis*](https://doi.org/10.1038/s41592-025-02909-7), Nature Methods 23, 350–359 (2026) | Pertpy 框架与实现引用。 |
| Milo 原论文 | Dann et al., [*Differential abundance testing on single-cell data using K-nearest neighbour graphs*](https://doi.org/10.1038/s41587-021-01033-z), Nature Biotechnology 40, 245–253 (2022); [开放全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7617075/) | 重叠 KNN 邻域、按实验 Sample 计数、负二项 GLM、edgeR QL F 检验、密度加权空间 FDR。 |
| miloR 官方实现 | [`testNhoods`](https://marionilab.github.io/miloR/reference/testNhoods.html)、[`graphSpatialFDR`](https://marionilab.github.io/miloR/reference/graphSpatialFDR.html)、[官方教程](https://marionilab.github.io/miloR/articles/milo_demo.html) | 原始 R 接口的 design/contrast、QL 检验、`k-distance` 空间 FDR、结果列和三对三 Sample 示例。 |
| edgeR/limma | [edgeR 包页](https://bioconductor.org/packages/release/bioc/html/edgeR.html)、[edgeR User's Guide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)、[edgeR manual](https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf)、[limma manual](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf) | NB GLM、TMM、QL 拟合、设计矩阵秩/残余自由度以及 contrast 是设计系数线性组合。 |
| Scanpy 邻居图 | [`scanpy.pp.neighbors`](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.neighbors.html) | `use_rep`、`n_neighbors`、`random_state`、`key_added` 的含义及 `.uns/.obsp` 存储约定。 |

配套统计引用应至少包括：Robinson, McCarthy & Smyth, [edgeR](https://doi.org/10.1093/bioinformatics/btp616)；Robinson & Oshlack, [TMM](https://doi.org/10.1186/gb-2010-11-3-r25)；Lun, Chen & Smyth, [QL edgeR workflow](https://doi.org/10.1007/978-1-4939-3578-9_19)。这些引用属于方法披露，不代替对本次设计、复制和混杂的说明。

## 2. Pertpy 1.3.0 精确接口

以 `scverse/pertpy` 的 `1.3.0` 标签源码为准，相关公开方法签名为：

```python
Milo.load(input, feature_key="rna")

Milo.make_nhoods(
    data,
    *,
    neighbors_key=None,
    feature_key="rna",
    prop=0.1,
    seed=0,
    copy=False,
)

Milo.count_nhoods(data, sample_col, feature_key="rna")

Milo.da_nhoods(
    mdata,
    *,
    design,
    model_contrasts=None,
    subset_samples=None,
    add_intercept=True,
    feature_key="rna",
    reml=True,
    max_iter=50,
    tol=1e-5,
    solver="pydeseq2",
)

Milo.annotate_nhoods(mdata, anno_col, feature_key="rna")
```

必须注意以下版本事实：

- `da_nhoods` 的 1.3.0 默认求解器是 `pydeseq2`，不是原始 Milo 的 edgeR QL 路径。官方文档称 `edger` 最接近 R 实现。因此科研节点不能省略 `solver` 并把结果统称为 edgeR/Milo QL。
- `subset_samples` 已弃用。官方源码要求需要子集时先对子细胞级数据，再重建邻居图和邻域；不能在已由全部细胞定义的图上事后删 Sample。
- 1.3.0 明确检查设计必须保留残余自由度，即 `n_samples > rank(design)`。这是“可估计”的必要条件，不是生物重复充分、无混杂或统计功效充分的保证。
- 随机截距会切换到另一套负二项混合模型；其估计量和结果列不同于 edgeR QL。重复测量不是在同一原子接口里添加一个字符串即可安全支持的细节。

## 3. 官方工作流及对象语义

### 3.1 `load`

`milo.load(adata, feature_key="rna")` 把细胞级 `AnnData` 放入 `MuData["rna"]`，并创建 Milo 结果模态。它是数据容器转换，不进行邻居计算或差异检验。

### 3.2 邻居图和表示空间

Pertpy 教程先在细胞模态调用 Scanpy，例如：

```python
sc.pp.neighbors(mdata["rna"], use_rep="X_scVI", n_neighbors=150)
```

`scanpy.pp.neighbors` 的 `use_rep` 指定用于距离计算的矩阵，`n_neighbors` 控制局部图尺度。若设置 `key_added="name"`，Scanpy 写入：

- `.uns["name"]`，含图参数和表示空间元数据；
- `.obsp["name_distances"]`；
- `.obsp["name_connectivities"]`。

Pertpy `make_nhoods(neighbors_key="name")` 正是按这套约定读取图；`neighbors_key=None` 读取默认 `neighbors/connectivities/distances`。Pertpy 还从 `.uns[neighbors_key]["params"]["use_rep"]` 取用于邻域精炼的表示；元数据缺失时可能回退至 `X_pca` 并警告。因此 `neighbors_key` 与 `representation` 不是两个可任意组合的独立科学参数，适配器必须构建并验证同一个私有图工件。

Milo 原论文假设 KNN 图忠实表达细胞表型流形；技术效应应尽量在构图前处理，也可在 Sample 级 GLM 中纳入剩余技术协变量。换言之，批次校正表示并不自动消除设计矩阵中 Sample 级批次项的必要性，两者解决的问题也不完全相同。

### 3.3 `make_nhoods`

`make_nhoods` 从 KNN 图随机抽取约 `round(n_obs * prop)` 个候选顶点，再把候选精炼至其邻域在指定低维表示中更具代表性的实际细胞，并去除精炼后重复的索引细胞。1.3.0 源码使用显式 `seed` 控制随机抽样，并写入：

- `obsm["nhoods"]`：cell × neighborhood 的二元稀疏成员矩阵；
- `obs["nhood_ixs_random"]` 和 `obs["nhood_ixs_refined"]`；
- `obs["nhood_kth_distance"]`；
- `uns["nhood_neighbors_key"]`。

`prop` 是候选索引顶点比例，不是最终邻域覆盖率或固定的最终邻域数；精炼去重会使最终数量小于候选数量。官方教程明确提示 `k` 和 `prop` 应按数据规模、供体数和邻域大小分布做敏感性考虑，而不是把某个教程值视为通用常数。

源码通过 Python `random.seed(seed)` 与 `random.sample` 实现抽样。这会改变进程级 Python RNG 状态；宿主节点若要求无全局副作用，必须在受控临界区保存并恢复 RNG 状态。

### 3.4 `count_nhoods`：统计单位是 Sample

`count_nhoods(..., sample_col=...)` 把每个邻域中的细胞按实验 Sample 计数，产生 Sample × neighborhood 计数矩阵。每个 Sample 是 GLM 的一行；细胞不是独立重复。官方教程和原论文都把 Sample/供体/生物学重复作为变异估计单位。

这要求：

- 每个 `sample_col` 值代表独立生物学标本，而不是 Condition 标签、单个细胞或技术孔；
- Condition 及每个协变量在一个 Sample 内必须唯一；
- 不能以细胞数大为由弥补 Sample 数不足；
- 同一供体的重复测量不是独立 Sample，需专门的重复测量模型，不能在独立 Sample 接口中假装独立。

### 3.5 `da_nhoods`：设计、contrast 与求解器

设计公式描述 Sample 级模型；contrast 是设计矩阵系数的线性组合。miloR `testNhoods` 在没有显式 contrast 时默认检验公式最后一项/设计矩阵最后一列，但多水平因子、协变量顺序和编码会使这个隐式规则不适合作为工作流节点契约。官方 edgeR/limma 文档同样要求设计矩阵满秩、保留残余自由度，并明确检验系数或数值 contrast。

Pertpy 1.3.0 的两个非混合求解器不是同一个统计方法：

- `solver="edger"`：TMM 标准化，edgeR 负二项 GLM，离散度估计，`glmQLFit(..., robust=True)` 与 `glmQLFTest`；复杂 contrast 经 limma `makeContrasts`。这是最接近 Milo 原始 R 实现的路径。
- `solver="pydeseq2"`：PyDESeq2 负二项 Wald 检验，并把结果重命名到部分 Milo 兼容列。它不是 edgeR QL F 检验，即使同一高层方法名下结果通常相近。

因此可复现接口必须显式固定求解器和 contrast 方向。若目标是复刻原始 Milo，应传 `solver="edger"`，并记录 R、edgeR、limma、statmod、rpy2 与 Pertpy 版本。

复制数方面，Pertpy 只硬性要求实际设计有正残余自由度；miloR 官方示例使用每个 Condition 3 个 Sample，原论文的平衡模拟也以 3 个重复/Condition 作为最低模拟设置并显示增加重复改善效应估计。OpenBio 若采用“每个比较 Condition 至少 3 个独立 Sample”，应明确标注这是保守的产品前置条件，而不是官方库宣称的普适功效阈值；满秩、正残余自由度、无完全混杂仍需单独检查。

### 3.6 空间 FDR

邻域重叠使普通 BH 假设不适合直接作为 Milo 的主要发现控制。Milo 原论文与 miloR `graphSpatialFDR` 使用密度/图结构加权的 BH；Pertpy 1.3.0 的 `SpatialFDR` 以 `1 / kth_distance` 为权重执行加权校正。`FDR` 是普通 BH，`SpatialFDR` 才应作为 Milo 报告中的主要多重检验量。

这也意味着 `kth_distance` 必须与本次图及索引细胞一一对应。零、非有限或错轴距离不能被静默修补成可信空间校正。1.3.0 对部分无穷权重/缺失校正结果有防御性填充值；严格科研适配器仍应在调用前后验证距离和结果，而不是把填充视为有效证据。

### 3.7 `annotate_nhoods`

`annotate_nhoods(..., anno_col=...)` 对每个邻域取细胞注释的多数类别，并写入：

- `var["nhood_annotation"]`：多数标签；
- `var["nhood_annotation_frac"]`：多数标签所占比例；
- `varm["frac_annotation"]`：各标签比例矩阵；
- `uns["annotation_labels"]` 与相关元数据。

源码拒绝数值注释，并通过 one-hot 编码计算比例。缺失标签不会作为一个显式类别进入分母，可能抬高表面纯度。因此适配器应在调用前拒绝缺失/空白注释；“低于阈值称 Mixed”应是额外派生的展示字段，不能覆盖官方多数标签和原始比例。

## 4. 结果列与解释

Pertpy edgeR 路径在 Milo 结果 `.var` 中产生或保留的核心列为：

| Pertpy 列 | 含义 | 稳定输出建议 |
| --- | --- | --- |
| `index_cell` | 邻域索引细胞身份/位置 | 保留明确的 observation ID，不以隐式行号代替。 |
| `kth_distance` | 索引细胞到第 k 邻居的距离，也是空间 FDR 权重依据 | 要求有限且大于 0。 |
| `logFC` | 模型比较的 log2 abundance fold change | 明确固定为 comparison 相对 reference；正值表示 comparison 富集。 |
| `logCPM` | 邻域计数的平均 log counts-per-million 尺度 | 披露但不当作效应方向。 |
| `F` | edgeR 准似然 F 统计量 | 仅 edgeR QL 固定路径输出。 |
| `PValue` | 未校正 QL 检验 p 值 | 披露完整。 |
| `FDR` | 普通 BH 校正 | 保留作诊断，不作为 Milo 的主要发现阈值。 |
| `SpatialFDR` | 按邻域密度加权的空间 FDR | 主要多重检验量。 |
| `nhood_annotation` | 多数细胞注释 | 原值保留。 |
| `nhood_annotation_frac` | 多数注释比例 | 原值保留，再派生 `is_mixed`。 |

邻域彼此重叠，一个细胞可属于多个邻域。因此“显著邻域数”不是独立细胞群数量，“注释为某 cell type”也不等于该 cell type 整体显著。报告必须把发现表述为图上局部表型区域的差异丰度证据。

## 5. 当前仓库依赖事实

改前检查发现仓库元数据锁定 Scanpy/AnnData，但没有声明 Pertpy、MuData、rpy2 或 R/Bioconductor Milo 运行栈；也没有 Milo 专项例程测试。Pertpy 官方安装资料表明 edgeR 路径不仅需要 Python 侧 `rpy2`，还需要外部 R 及 edgeR/limma/statmod。后续实现不能把 import 成功等同于完整后端可用，必须给出可操作的依赖错误并记录实际后端版本。

## 6. 官方用法对节点的约束结论

1. 执行顺序应是：验证并按明确的两 Condition 选择细胞 → 在选中细胞副本上构图 → 随机且可复现地抽样/精炼邻域 → 按独立 Sample 计数 → 构建设计与显式 comparison-reference contrast → 显式 `solver="edger"` 运行 QL → 空间 FDR → 描述性注释 → 规范化结果与报告。
2. `Sample` 是统计重复单位；每个 Condition 至少 3 个独立 Sample 是 OpenBio 的保守最低门槛，同时仍需满秩、正残余自由度与协变量不完全混杂。
3. 表示空间、图、索引细胞、邻域成员、Sample 计数和 `kth_distance` 是一个不可拆开的证据链；不能接受陈旧或外来图键而不验证来源。
4. 不能省略求解器，也不能把 PyDESeq2 Wald 结果标成 edgeR QL。
5. `SpatialFDR` 是主要校正量；普通 `FDR` 只作补充。
6. 注释是对已检验邻域的描述性解释，不参与 DA 模型；原多数标签/比例必须保留，provisional/curated 状态必须披露。
7. 一般科研报告至少应披露 Sample 数与每组复制、设计矩阵和 contrast 方向、表示/图参数、随机种子、邻域数及大小、求解器和软件版本、全部检验数、空间 FDR 阈值、显著方向、注释纯度以及方法局限。
