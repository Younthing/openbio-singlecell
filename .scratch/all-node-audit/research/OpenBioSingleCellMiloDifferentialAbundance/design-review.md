# OpenBioSingleCellMiloDifferentialAbundance：深模块设计评审

本文是代码修改前的设计记录。结论以同目录 `official-usage.md` 核对的 Pertpy 1.3.0、miloR、Milo 论文与 edgeR 官方资料为依据。

## 1. 当前节点与集成现状

当前 `OpenBioSingleCellMiloDifferentialAbundance` 接收：

```text
adata
sample_key="sample"
design="~ group"
contrast=""
annotation_key="cell_type"
use_rep="X_pca"
n_neighbors=30
neighborhood_proportion=0.1
mixed_annotation_threshold=0.6
neighbors_key="openbio_milo"
random_seed=123
```

它复制 AnnData，依次执行 `Milo.load`、`scanpy.pp.neighbors`、`make_nhoods`、`count_nhoods`、`da_nhoods` 和 `annotate_nhoods`，再返回原始 Milo `.var`。当前只有 `table` 输出。

仓库检查同时发现：

- 注册测试只证明类可注册，没有 Milo 科学行为专项测试；
- 当前示例工作流没有连接该节点；
- 项目依赖没有声明 Pertpy/MuData/rpy2/R/edgeR 运行栈；
- 没有可证明旧工作流究竟使用 edgeR 还是 Pertpy 默认 PyDESeq2 的集成工件。

因此现状不是“接口已稳定但报告不足”，而是统计 estimand、后端与可复现性都未冻结。

## 2. 正确性审计

### P0：模型重复单位没有被保护

节点只检查 `sample_key` 列存在。它不能阻止用户把 Condition、技术孔或单细胞 ID 当作 Sample，也不检查一个 Sample 内 Condition/协变量是否唯一。Milo 的模型行是按 Sample 聚合的邻域细胞数；把细胞或 Condition 当重复会分别造成伪重复或完全没有方差估计。

### P0：自由文本 design/contrast 不是安全的原子 estimand

`design="~ group"` 与自由文本 `contrast` 暴露了上游公式编码、因子参考、R 合法列名和系数顺序，却没有输出设计列、秩、残余自由度或 contrast 向量。多水平因子、缺失水平、交互项、变量名转义、条件与批次完全混杂都可能得到错误方向或不可估计模型。

一般节点应只表达一个清楚的问题：在选定 reference 与 comparison Condition 间，调整声明的 Sample 级加性 nuisance covariates 后，哪些图邻域差异丰度？交互、连续暴露、多重 contrast 和重复测量应由另一个明确命名的高级设计模块处理。

### P0：后端默认值已发生统计语义漂移

当前调用没有传 `solver`。Pertpy 1.3.0 默认 `pydeseq2`，而原始 Milo 使用 edgeR QL F 检验；两者不是相同检验。节点名称和当前实现无法证明历史结果属于哪条路径，也没有报告后端。新节点必须显式固定原始 Milo 路径 `solver="edger"`，不能依赖第三方默认值。

### P0：`random_seed` 没有传给真正随机的步骤

当前 seed 只传给 `scanpy.pp.neighbors`，但 `Milo.make_nhoods` 的随机采样有独立 `seed` 参数，当前未传，实际使用 Pertpy 默认 0。UI 中的 123 因此不能完整描述结果。源码还会重置进程级 Python `random` 状态；并发宿主中这是额外副作用。

### P1：图来源可被误配，公开 `neighbors_key` 增加无效状态

当前节点自己构图，因此用户提供 `neighbors_key` 只是内部存储名，不是科学自由度。公开它会引入与现有 `.uns/.obsp` 碰撞、错读旧图和表示元数据不一致的风险。图应在选定两个 Condition 的私有副本上构建，键由适配器生成和验证，用户只选择科学上有意义的 `representation_key`、`n_neighbors` 与邻域比例。

### P1：未在构图前限定比较群体

若输入含第三个 Condition，当前图和邻域由所有细胞共同定义，而模型 contrast 可能只讨论其中两个水平。Pertpy 官方已弃用事后 `subset_samples`，要求先子集细胞再重建图。成对节点应在任何图操作前选择 reference/comparison 两组，并报告被排除的 Condition/细胞数。

### P1：注释缺失与 Mixed 覆盖破坏证据

当前节点不验证注释缺失/空白/数值类型。Pertpy 的 one-hot 计算会从比例分母隐式忽略缺失标签。节点随后把低纯度的 `nhood_annotation` 原地改成 `Mixed`，丢失官方计算出的多数标签。正确输出应同时保留原多数标签、比例，并额外给出 `reported_annotation` 与 `is_mixed`。

### P1：没有图、计数、设计和结果后置条件

当前没有检查：表示有限性、observation 轴、KNN 图键、成员矩阵二元性、索引细胞唯一性、正 `kth_distance`、每个 Sample 计数行、整数非负计数、非零 library、设计满秩、contrast 可估计、结果一行一邻域、数值范围和 SpatialFDR 完整性。第三方版本漂移或畸形输出可能被直接交付。

### P1：输出既不稳定也不足以形成报告

原始 `.var` 列随后端而变，列名没有明确 comparison/reference 方向；普通 `FDR` 与 `SpatialFDR` 没有区分；没有 Sample 复制、图随机性、测试总数、显著方向、注释状态、引用、软件版本、严格 JSON `summary` 或等效 `code`。

## 3. 模块边界备选方案

### 方案 A：把 load、构图、抽邻域、计数、检验和注释拆成公开节点

优点是每一步看似可复用；代价是把 MuData 槽位、私有键、轴一致性、seed 和步骤顺序全部暴露给工作流。用户可将全体 Condition 的旧图与新子集计数拼接，或跳过必要阶段。删除任一薄包装节点也不会损失领域抽象，只会少一个调用转发，说明这些不是有深度的公开模块。

结论：不采用。阶段可以是私有函数，不能成为无契约的公共中间节点。

### 方案 B：保留一个“任意 design/contrast 的通用 Milo”节点

优点是表面灵活；代价是接口泄露公式语言、R 命名和求解器实现，节点不能知道用户要检验的 estimand，也无法稳定生成报告或迁移。它对一般用户是浅模块，对高级用户又不够完整，因为没有显式设计工件和 contrast 矩阵。

结论：不采用。未来若确需高级设计，应输入一个经过独立验证、可审计的 Sample-design artifact，而不是任意字符串。

### 方案 C：一个成对 Sample-level Condition Milo 深模块

模块内部拥有：精确的两组细胞选择、私有图、可复现邻域、Sample 计数、安全设计、一个 comparison-reference contrast、edgeR QL、空间 FDR 与邻域描述。公开接口只暴露改变科学问题或合理敏感性分析的参数。

结论：采用。它以较小接口隐藏一条较复杂但内聚的证据链，具有更高 depth/leverage 和更好的 locality。

## 4. 保留、合并与删除决定

**保留并破坏性重构节点概念**，其原子操作定义为：

> 在一个明确的细胞表示上，为恰好两个 Condition 的细胞构建 Milo 邻域；以独立 Sample 为行、调整声明的加性 nuisance covariates，运行一个固定方向的 edgeR QL abundance contrast；对全部邻域执行空间 FDR，并附加描述性注释及完整报告。

不与 cluster composition、scCODA 或 cluster-level pseudobulk 合并。Milo 的连续图邻域、重叠假设和空间 FDR 是独立方法语义；合并会把“离散群比例”和“图上局部丰度”混为一谈。

不把 `annotate_nhoods` 单独做成公共节点。它依赖本次成员矩阵且只为解释本次结果服务，留在模块内部更内聚；原始多数标签和派生 Mixed 状态都输出即可。

可与 pseudobulk Condition inference 节点共享一个**私有、纯函数**的 Sample-role/design preflight seam：角色解析、Sample 内唯一映射、编码、秩、残余自由度、confounding 与 contrast estimability。该 helper 不应知道 Milo 图或 edgeR 结果列。Milo 自己拥有构图、邻域、计数、空间 FDR 和注释适配器。

删除下列旧接口概念：自由文本 `design`、自由文本 `contrast`、公开 `neighbors_key`、后端默认继承、覆盖原始注释、只返回上游 `.var`。删除它们会减少无效状态，不会损失目标科学能力。

## 5. 拟议原子接口

### 输入顺序与默认值

```text
adata: OPENBIO_ANNDATA
sample_key: STRING = "sample"
condition_key: STRING = "condition"
reference_condition: STRING                  # 必须显式非空
comparison_condition: STRING                 # 必须显式非空
technical_batch_key: STRING = ""            # 可选 Sample 级分类 nuisance
categorical_covariate_keys_json: STRING = "[]"
continuous_covariate_keys_json: STRING = "[]"
annotation_key: STRING = "cell_type"
annotation_status: COMBO["provisional", "curated"] = "provisional"
representation_key: STRING = "X_pca"
n_neighbors: INT = 30
neighborhood_proportion: FLOAT = 0.1
mixed_annotation_threshold: FLOAT = 0.6
spatial_fdr_threshold: FLOAT = 0.1
min_abs_log2_fold_change: FLOAT = 0.0
random_seed: INT = 123
```

两个 `*_keys_json` 必须是严格 JSON 字符串数组：元素仅允许 canonical 非空字符串、不得重复，且不得与 Sample/Condition/Technical batch/annotation 角色重叠。使用 JSON 而不是逗号拼接可避免列名中的逗号与空白歧义。分类 nuisance 参考水平按原 categorical 声明顺序；非 categorical 列按所选 Sample 的首次出现顺序，实际顺序必须写入 summary。连续协变量保持原始单位，记录范围，不做隐式 winsorization 或缺失插补。

### 应暴露的参数

- Sample、Condition、reference/comparison：定义 estimand 和独立重复；必须可见。
- Technical batch 与有限的加性分类/连续 covariates：改变设计；作为高级参数可见。
- representation、`n_neighbors`、`neighborhood_proportion`、seed：改变被检验的邻域；必须可见并进入方法摘要。
- annotation 及其 provisional/curated 状态、Mixed 阈值：改变解释标签，不改变模型；必须区分报告层与检验层。
- SpatialFDR/LFC 阈值：只改变 calls 与摘要，不过滤 table，也不改变拟合/校正全集。

### 应隐藏且固定的参数

- `feature_key="rna"`、私有 collision-free `neighbors_key`；
- 先按两个 Condition 子集，再构图；
- Scanpy 图策略固定并披露为 `metric="euclidean"`、`knn=True`、`method="umap"`、`transformer="pynndescent"` 与 `random_state=random_seed`；不使用会按数据规模改变算法的 `transformer=None` 自动选择；
- Pertpy `make_nhoods(..., prop=..., seed=random_seed, copy=False)`；
- `count_nhoods(..., sample_col=<private safe alias>, feature_key="rna")`；
- 安全内部列名和安全水平 `reference`/`comparison`，无截距、Condition 第一项的加性设计，显式 comparison-minus-reference contrast；
- `da_nhoods(..., add_intercept=False, solver="edger", feature_key="rna")`，不使用弃用的 `subset_samples`；
- edgeR TMM、robust QL 拟合/QLF 检验和 `k-distance` 空间 FDR；
- 全部邻域进入多重检验，不按 p 值/top-N 截断；
- 不修改输入，不自动回退 PyDESeq2、混合模型、cluster-level 检验或 cell-level 检验。

若 Pertpy 公开接口不能稳定表达固定策略，应在一个窄的带版本断言 adapter 中调用受支持接口，或清楚失败；不得假称执行 edgeR QL 而实际继承别的默认值。

### 输出顺序

```text
table: OPENBIO_TABLE
summary: STRING       # 严格 JSON object
code: STRING          # 科学等效、可执行 Python 源码
```

该节点不输出修改后的 AnnData/MuData，因为它的科学产物是一次 contrast 的邻域结果；暴露内部 MuData 会扩大耦合并允许下游误用未冻结的 Pertpy 槽位。

## 6. 统计状态与允许的声明

满足全部前置条件时，本节点可标记为 `sample_level_condition_inference`：模型行是独立 Sample，离散度来自 Sample 间变化，且使用预先声明的 comparison-reference contrast。它不是 cell-level inference，也不因细胞数量多而增加 n。

最低产品门槛固定为每个所选 Condition 至少 3 个独立 Sample。这比 Pertpy 的“正残余自由度”检查更严格，是基于原论文/官方示例的保守最低实践门槛，而非功效保证。仍必须要求设计满秩、正残余 df、contrast 可估计且不存在完全混杂。

若输入是配对/纵向/同一供体重复 Sample，当前原子接口必须拒绝，并指向未来重复测量 Milo 模块；不得把重复记录当独立 Sample。若 annotation 是 `provisional`，DA 推断仍可成立，但标签解释必须始终称 provisional。无随机分配或其他因果设计证据时，摘要使用“associated with / 相较于”而不是“caused by / 导致”。

## 7. 输入验证与原子性

### 7.1 身份与角色

- `obs_names` 唯一、非空、无前后空白；所有角色列存在且轴长等于 `n_obs`。
- Sample、Condition、分类 nuisance 和 annotation 只允许非缺失 canonical 字符串；拒绝 `1` 与 `"1"` 之类的强制转字符串碰撞。
- reference/comparison 均存在、互异且精确匹配；记录并排除其他 Condition 的细胞。
- 每个 Sample 映射到恰好一个 Condition、一个 Technical batch 和每个分类/连续协变量值；一个 Sample 内出现多个值立即失败。
- Sample、Condition、batch、covariate、annotation 列角色互不别名；reference/comparison 各至少 3 个独立 Sample。
- 连续 covariates 必须在 Sample 级有限、非布尔、非恒定；分类 covariates 至少两个实际水平，但不能一 Sample 一水平。

### 7.2 表示与图参数

- `representation_key` 存在于 `.obsm`，二维，行数精确等于 `n_obs`，至少 2 列，数值、有限且不是所有行相同；记录 shape、dtype、内容 fingerprint 和已知生成元数据。
- 本方法不读取表达矩阵作为计数输入；它使用 representation 定义表型图、使用细胞归属形成 Sample 邻域计数。因此不得伪称验证了 raw counts/CP10K/log1p 表达状态，必须披露“representation-derived graph”。
- `2 <= n_neighbors < selected_n_cells`；`0 < neighborhood_proportion <= 1` 且 `round(n_cells * prop) >= 1`；seed 是支持范围内整数；阈值有限并位于合法区间。
- 选中细胞数、每 Sample 细胞数和每 Condition 细胞数全部披露。极端不平衡发出警告，但不能用细胞数代替 Sample 门槛。

### 7.3 设计矩阵

- 在安全内部别名上确定性编码 additive design；Condition 固定两水平并显式生成 comparison-reference contrast。
- 设计行顺序必须与 Sample 计数矩阵一致；列名、分类参考水平、数值范围、矩阵 rank 和 `residual_df = n_samples - rank` 写入 summary。
- 要求 rank 等于列数、residual df > 0、contrast 非零/有限且位于设计可估计空间；报告具体冲突项后拒绝完全混杂，不自动删协变量。
- 仅包含选中两 Condition 的 Sample。任何后端静默删除零 library Sample 都视为契约失败；应在拟合前发现。

### 7.4 每阶段后置条件

- Scanpy 图的 `.uns/.obsp` 私有键完整、shape 为 cell × cell、轴与选中 observation 身份一致、距离/连接有限非负；图元数据的 representation 和 k 与请求一致。
- 邻域成员矩阵为 cell × nhood 二元矩阵；至少产生两个唯一邻域；随机/精炼索引、成员列、索引细胞和 `kth_distance` 一一对应；每个邻域非空，距离有限且严格大于 0。
- Sample × nhood 计数为非负整数；包含每个且仅包含选中 Sample；每个 Sample library total 大于 0；列与成员矩阵邻域完全一致。
- annotation 调用前无缺失/空白/数值标签；官方多数标签、比例轴和范围完整。
- edgeR 结果恰好一行/邻域，无重复/丢失/额外行；`logFC/logCPM/F/PValue/FDR/SpatialFDR` 均有限，`F >= 0`，所有概率在 `[0, 1]`。固定 edgeR QL 路径若不能给出完整合法结果，整次执行失败并说明邻域，而不是静默丢行或用第三方填充值冒充证据。

### 7.5 失败原子性与全局状态

所有第三方调用只作用于私有深拷贝；输入的 observation/variable identity、role 列和 representation fingerprint 在异常和成功后都不得改变。完成全部后置条件前不发布任何输出。

由于 Pertpy 1.3.0 邻域抽样会调用进程级 `random.seed`，实现需在进程锁内保存/恢复 Python `random` 状态，并对实际使用的 NumPy/Scanpy RNG 采用显式种子或同等隔离；`finally` 中恢复，即使后端抛错也不能污染其他节点。不得改写全局 R options 或默认 contrasts；如桥接层需要临时状态，同样做受控恢复。

## 8. canonical `table` 契约

结果按确定的邻域身份顺序输出完整全集，不按显著性排序或截断：

```text
neighborhood_id: string
index_cell: string
neighborhood_size: integer
kth_distance: float
reference_condition: string
comparison_condition: string
log2_fold_change: float
log_counts_per_million: float
quasi_likelihood_f: float
p_value: float
p_adjusted_bh: float
spatial_fdr: float
majority_annotation: string
majority_annotation_fraction: float
reported_annotation: string
is_mixed: boolean
```

`neighborhood_id` 是稳定、唯一、无歧义的适配器 ID；`index_cell` 保留原 observation ID。`log2_fold_change > 0` 始终表示 comparison 相较 reference 更丰富。`reported_annotation` 在多数比例低于 `mixed_annotation_threshold` 时为 `Mixed`，否则等于 `majority_annotation`；原标签和比例永不覆盖。阈值只用于 calls/summary，不删除 table 行，也不改变 SpatialFDR 的假设全集。

普通 BH 与空间 FDR 必须同时保留并明确区分。主结论只用 `spatial_fdr <= spatial_fdr_threshold`，再按可选 `abs(log2_fold_change) >= min_abs_log2_fold_change` 形成报告 calls。

## 9. `summary` JSON 契约

`summary` 必须是标准 JSON object 文本，使用 `allow_nan=False`；不得含 Python repr、NaN 或 Infinity。建议 schema：

```json
{
  "schema_version": "openbio-singlecell/milo-differential-abundance/v1",
  "node_id": "OpenBioSingleCellMiloDifferentialAbundance",
  "analysis_status": "sample_level_condition_inference",
  "method": {},
  "comparison": {},
  "input_evidence": {},
  "design_evidence": {},
  "graph_evidence": {},
  "analysis_summary": {},
  "key_results": [],
  "parameters": {},
  "warnings": [],
  "limitations": [],
  "references": [],
  "software_versions": {}
}
```

最低披露内容：

- 方法：Milo overlapping KNN neighborhoods、Sample-level counts、TMM、edgeR robust QL F、普通 BH 与 k-distance SpatialFDR 的清楚区别；
- comparison：Condition key、原始 reference/comparison 标签、正 LFC 方向；
- 输入证据：输入/选中/排除细胞数，Sample ID/Condition 映射摘要，每组独立 Sample 数，每 Sample/Condition 细胞数，annotation 状态，representation key/shape/fingerprint/已知 provenance；
- 设计证据：固定公式的可读渲染、实际设计列与行顺序、分类参考、rank、residual df、contrast 向量、confounding/estimability 检查；
- 图证据：Scanpy 参数、邻域 seed、请求/实际抽样数、最终邻域数、覆盖率、邻域大小和 kth-distance 分布、每邻域 Sample/Condition 覆盖诊断；
- 分析摘要：总检验数、合法数、SpatialFDR 阈值下 comparison-enriched/reference-enriched 数、Mixed 数、无显著结果时的明确空结论；
- `key_results`：有界的最强上调/下调显著邻域记录，至少含 ID、index cell、LFC、SpatialFDR、多数标签/比例与大小；它们来自完整 table，不构成另一次筛选检验；
- 版本：OpenBio、Python、Pertpy、Scanpy、AnnData、MuData、NumPy、pandas、SciPy、scikit-learn、PyNNDescent、rpy2、R、edgeR、limma、statmod；缺失版本视为不可完整披露；
- 引用：Milo 原论文、Pertpy 论文、edgeR、TMM、QL workflow，使用结构化 list；
- 限制：结果依赖 representation/k/prop/seed；重叠显著邻域不是独立 cell type；标签只是描述；批次校正表示可能不完整；三 Sample/组只是最低门槛而非功效保证；观察研究不能作因果解释；细胞捕获深度与组成性会影响 abundance。

报告语言必须区分 Sample n 与 cell count，区分普通 BH 与 SpatialFDR，区分 provisional annotation 与 curated annotation。没有显著邻域时必须直接报告 0，而不能把最小 p 值写成阳性发现。

## 10. `code` 科学等效契约

`code` 输出一段可执行 Python 源码，定义例如：

```python
def run_milo_differential_abundance(adata, ...):
    ...
    return table, summary_dict
```

它必须包含与运行节点相同的：

- strict JSON 参数解析和全部身份/角色/表示验证；
- 先选择两个 Condition、再在副本构图；
- 相同 Scanpy 图参数和相同 `make_nhoods(seed=...)`；
- RNG 状态隔离；
- Sample 计数、safe alias、设计矩阵/秩/contrast 检查；
- 显式 `solver="edger"`、`add_intercept=False` 及其余固定公开参数；
- annotation 原值保留、canonical 列映射、完整后置条件；
- 相同 calls、summary 字段、引用和软件版本读取。

源码不得依赖 ComfyUI 节点对象、隐藏文件、网络下载、预存 MuData 或未披露全局变量。它可以省略 UI timing/插件历史，但不能省略会改变科学结果、拒绝条件或报告解释的逻辑。测试应在同一输入/backend 下比较运行节点与 `exec(code)` 的完整 table、设计/contrast、邻域身份和关键 summary 字段，而不仅是能编译。

## 11. 依赖、版本漂移与迁移

### 依赖策略

后续实现需要显式的可选 Milo extra，并在运行前按顺序探测 Pertpy、MuData、rpy2、R runtime、edgeR、limma、statmod。错误要给出检测到的版本和安装动作。Pertpy 的 Python extra 不能代替外部 R/Bioconductor 安装。不得静默回退 PyDESeq2，因为这会改变方法。

建议冻结并运行时核验一个经过测试的窄版本窗口，例如 Pertpy `>=1.3,<1.4`，同时检查上述五个公开 Milo 方法签名。版本不兼容时失败，而不是按反射猜参数。变更记录显示 1.0 引入 PyDESeq2、1.2 引入混合模型、1.3 增加残余自由度错误；未来 Scanpy 文档也可能改变 RNG 命名，证明浮动默认不可接受。

### 工作流迁移风险

- 输出从一个 `table` 增至 `table, summary, code`；旧连线索引必须迁移。
- `design`、`contrast`、`neighbors_key` 被删除，`use_rep` 改名为 `representation_key`；新增明确 Condition/reference/comparison 和 annotation 状态。
- 旧 `random_seed` 实际未控制 Milo 抽样；迁移后同值不保证复现旧邻域。
- 历史节点因 Pertpy 版本/默认值可能实际运行 PyDESeq2，也可能在另一环境运行不同路径。没有旧 summary 就无法可靠推断；迁移不能伪装为数值兼容，应提示用户重新审核和运行。
- 新节点要求每组至少 3 个独立 Sample，并拒绝复杂公式、重复测量、完全混杂和第三 Condition 参与图；部分旧工作流将被有意阻止。
- 新的私有两组图与旧的全数据图定义不同；结果变化是修正后的 estimand，不是应隐藏的格式差异。
- 仓库目前无 Pertpy/edgeR 依赖声明和 Milo 示例；实现批次需同步注册、generator、迁移、依赖和示例，但本调研任务不修改这些文件。

推荐对旧节点 schema 使用明确的 migration-required shim 或有版本的迁移器；不要把任意公式自动翻译成 pairwise roles。只有能证明是简单两水平加性设计、能确定 reference/comparison 且能确定实际后端的实例才可机械迁移，其余要求人工确认。

## 12. 实现前测试计划

### 接口和输入

- 精确 schema、输入顺序/默认值、三输出端点及旧 schema migration-required 测试；
- 缺列、角色重叠、unknown Condition、同名 reference/comparison、缺失/空白/非字符串/强制转换碰撞；
- 一个 Sample 多 Condition/批次/协变量值，Condition 被错误用作 sample，少于 3 Sample/组；
- 分类/连续 covariate 的严格 JSON、重复 key、恒定/每 Sample 一水平、非有限数值、满秩/缺秩/零 residual df/不可估计 contrast；
- 配对或重复供体输入给出行动性拒绝，而不是当独立 Sample。

### 图和随机性

- representation 缺失、错轴、非二维、非数值、NaN/Infinity、退化；k/prop/seed 边界；
- fake Scanpy 验证只在两组副本构图、私有 key、representation、k、metric、transformer 和 seed 均显式；
- fake Pertpy 验证 `load → make_nhoods(seed=...) → count_nhoods → da_nhoods(solver="edger") → annotate_nhoods` 的顺序和全部科学参数；
- 两次同 seed 邻域身份一致，不同 seed 有可解释变化；成功和异常路径都恢复 Python RNG；并发执行不串扰；
- 畸形成员矩阵、重复索引、零/非有限 kth distance、错轴图、零 Sample library 均被拒绝。

### 设计、模型与方向

- 平衡三对三、带不混杂 batch、分类/连续 nuisance fixtures；Condition-batch 完全混杂、rank deficient 和 residual df 0 fixtures；
- 人工计数 fixture 证明正 `log2_fold_change` 是 comparison/reference；设计行、列、contrast 与 Sample 计数轴完全一致；
- fake backend 断言从不省略 solver、不使用 `subset_samples`、不调用 mixed/PyDESeq2 fallback；
- 可选的真实 edgeR 数值 smoke test 只在明确安装兼容栈时运行，并核对 TMM/QL/SpatialFDR 基本结果。

### 结果、报告与等效源码

- 缺列、重复/丢失/额外邻域、非有限统计量、负 F、越界 p/FDR、普通 BH 与 SpatialFDR 误映射均失败；
- missing/numeric annotation 失败；低纯度行保留多数标签并另加 Mixed；category order 和 observation identity 保持；
- table 永不按显著性截断，报告阈值不改变模型或多重检验全集；零显著/单向/双向结果摘要正确；
- summary 是 strict JSON，无 NaN/Infinity，包含完整 Sample/design/graph/方法/引用/版本/限制；provisional 语言始终保留；
- `code` 可编译、无 UI/网络依赖，并与节点在 canonical table、邻域、设计/contrast、calls 和 summary 核心字段上等价；
- 输入在成功、预检失败和 backend 失败后均逐槽不变，且不发布部分输出。

## 2026-08-29 report-contract repair

The node remains one atomic Sample-level contrast. Its method disclosures are now one report-ready `methods` string,
and all numeric/call highlights live in the canonical `key_results` mapping. The redundant `analysis_summary` branch
was removed; detailed comparison, design, graph, references, limitations, and versions remain additive report fields.
