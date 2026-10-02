# WP0 覆盖度核查报告（2026-10-01）

> 依据：《最终执行方案_二区落地_2026-10-01.md》WP0 条目。
> 结论一句话：**WP4 / WP5 / WP2 全部可行，无需缩减范围；仅 sapanisertib 为单库覆盖（PRISM），需按 AUC 为主口径预案执行。**

## 一、总判定表

| 核查项 | 结果 | 判定 |
|---|---|---|
| 5 个化合物在 GDSC1/GDSC2 的 IC50 覆盖 | 4/5 有（sapanisertib 两库均无） | ✅ 可行 |
| 5 个化合物在 PRISM 的覆盖 | 5/5 全有 | ✅ 可行 |
| DepMap 23Q4 Chronos 中 IGF2BP3 可用性 | 列存在，1100 模型全有数据 | ✅ 可行 |
| TCGA MC3 的 KEAP1/NFE2L2 注释 | 可用；KEAP1 突变 96/585（16.4%） | ✅ 可行 |
| LUAD 细胞系规模 | GDSC ~61、PRISM lung 102、Chronos 52 | 与方案预估 50–70 一致 |

## 二、5 个化合物的逐库覆盖（LUAD 口径）

| 化合物 | GDSC1（LUAD 有 IC50） | GDSC2 | PRISM 次级筛选（LUAD 模型 / 其中有 IC50） | 结论 |
|---|---|---|---|---|
| Torin-2 | 56 + 10（复测批次 ID 1202） | 无 | HTS002：44 / 44 | ✅ 双库 |
| BEZ235（Dactolisib / NVP-BEZ235） | 61 | 62 | HTS002：39 / 36；MTS010：42 / 13 | ✅ 三库 |
| MK-2206 | 56 | 62 | HTS002：40 / 2；MTS006：46 / 37 | ✅ 三库 |
| ZSTK474 | 61 | 无 | HTS002：43 / 33 | ✅ 双库 |
| Sapanisertib（MLN0128 / INK128） | 无 | 无 | HTS002：39 / 39；MTS010：45 / 28 | ⚠️ 仅 PRISM 单库 |

**执行建议（WP4 口径）**：
1. 逐药优先用覆盖更全的库做主分析，另一库做复现；sapanisertib 只报 PRISM，措辞按方案红线 2/6 处理。
2. PRISM 存在两个筛选批次：`MTS010` 为官方推荐批次（README 明确"recommended to use MTS010 when available"），但其 IC50 缺失明显多于 HTS002（如 BEZ235：MTS010 仅 13/42 有 IC50，HTS002 有 36/39）。**建议 IC50 与 AUC 双口径并列，或主用 AUC**（剂量反应曲线参数文件两者齐全），并在方法学中写明批次选择规则。
3. GDSC 主文件 `GDSC*_fitted_dose_response_24Jul22.csv` 来自 current_release（release 8.4 体系），LUAD 按 `TCGA_DESC == 'LUAD'` 过滤（该字段取值为 TCGA 缩写码，不是 "lung"）。

## 三、DepMap Chronos（WP5 输入核查）

- 文件：DepMap Public 23Q4 `CRISPRGeneEffect.csv`（Chronos 处理后基因效应，400,960,609 字节，与官方记录一致），1100 个细胞系模型。
- LUAD 口径：`Model.csv` 中 `DepmapModelType == 'LUAD'` 共 **89** 个模型；有 Chronos 数据的 **52** 个（与方案"LUAD 仅 50–70 个"的预估一致）。
- 关键基因 Chronos（LUAD 52 模型内）：

| 基因 | 中位数 | 最小值 | ≤ −1（强必需） | 判读 |
|---|---|---|---|---|
| IGF2BP3 | −0.061 | −0.344 | 0/52 | 非普遍必需，无单基因区分度 |
| SLC7A11 | +0.039 | −0.487 | 0/52 | 同上 |
| GCLC | −0.074 | −0.652 | 0/52 | 同上 |
| GCLM | +0.032 | −0.317 | 0/52 | 同上 |
| TFRC | −0.745 | −1.921 | 19/52 | 唯一接近泛必需者 |

**对修订 v2 的支撑**：数据印证了"不看单基因绝对必需性"的设计——除 TFRC 外四个基因在 LUAD 均无强依赖性，若按旧设计（单基因绝对值）必然全部阴性且无信息；分组比较（高/低防御组 CERES 差异，重点 IGF2BP3、SLC7A11）是唯一有区分度的口径。TFRC 在组间比较中预计也无差异（两家都接近地板效应），写结果时如实即可。

## 四、TCGA MC3 突变注释（WP2 输入核查）

来源：GDC API（TCGA 体细胞突变调用即 MC3 体系），TCGA-LUAD 共 585 例：

| 基因 | 突变病例数 | 占比 | 方案预估 | 核对 |
|---|---|---|---|---|
| KEAP1 | 96 / 585 | 16.4% | ~17% | ✅ 一致 |
| NFE2L2 | 17 / 585 | 2.9% | 未预估 | 样本偏少，亚组 KM 功效存疑 |
| TP53 | 281 / 585 | 48.0% | — | 供 WP7 突变景观用 |
| KRAS | 148 / 585 | 25.3% | — | 同上 |

**注意**：NFE2L2 突变仅 17 例，KEAP1/NRF2 联合亚组勉强可行（并组后 ~113 例）；单独做 NFE2L2 突变亚组 KM 功效不足，建议 WP2 中 KEAP1 与 NFE2L2 合并为"KEAP1–NRF2 通路突变"口径（豆包原建议也是合并写法），单基因口径放补充。

## 五、数据获取通道备忘（本机网络实测）

以下结论来自 2026-10-01 本机实测，后续 WP 复用：

| 通道 | 状态 | 说明 |
|---|---|---|
| depmap.org / cancerrxgene.org / cog.sanger.ac.uk / zenodo.org / api.figshare.com | ❌ 被 Cloudflare 或网络策略拦截 | 门户页、官网直链均不可达 |
| **ndownloader.figshare.com** | ✅ 直链可用 | DepMap/PRISM 官方 figshare 文件的下载通道；文件清单用 api.figshare.com 的 JSON（需经网页抓取工具代取） |
| **ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/** | ✅ 可用 | GDSC 官方 current_release，含 CSV 版 fitted dose-response |
| **api.gdc.cancer.gov** | ✅ 可用 | MC3 突变、病例计数；analysis 端点默认返回 TSV |
| www.orcestra.ca（API） | ✅ 可用 | 但其 PSet 下载链接指向 Zenodo（被拦），仅作目录查询 |

关键文件出处（写入方法学可引用）：
- GDSC1/2 fitted dose response：`ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC{1,2}_fitted_dose_response_24Jul22.csv`（release 8.4，2022-07-24）
- PRISM 次级筛选剂量反应参数：figshare 9393293 v4（Broad DepMap 官方，DOI 10.6084/m9.figshare.9393293.v4），`secondary-screen-dose-response-curve-parameters.csv`
- DepMap 23Q4：figshare 24667905，`CRISPRGeneEffect.csv`（Chronos）、`Model.csv`
- TCGA-LUAD 突变计数：GDC API（ssm_occurrences 端点，project TCGA-LUAD）

## 六、对后续 WP 的落地建议

1. **WP4（提前执行）**：按 §二 逐药定主库；sapanisertib 标注"探索性、单库"；IC50/AUC 双口径；LUAD n≈40–60/药，功效处于"探索性"档位，与方案措辞一致。
2. **WP5**：直接以已下载的 23Q4 `CRISPRGeneEffect.csv` + `Model.csv` 开工；分组比较 IGF2BP3、SLC7A11（另可对全 9 基因逐一比较，TFRC 预期阴性）。
3. **WP2**：KEAP1/NFE2L2 合并突变口径；MC3 注释若需样本级 MAF，可走 GDC 下载 `mc3.v0.2.8.PUBLIC.maf.gz` 或复用 cBioPortal luad_tcga_pan_can_atlas_2018（本次 API 503，下载时需重试）。
4. 已下载的 4 个原始文件保留在 `WP0_coverage/tmp/`（GDSC1/2、PRISM 次级筛选参数、CRISPRGeneEffect、Model、化合物表），WP4/5 直接复用，不必重复下载。
