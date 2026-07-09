# primerblast-oss

[![Release](https://img.shields.io/github/v/release/light-suzuki/primerblast-oss?sort=semver)](https://github.com/light-suzuki/primerblast-oss/releases)
[![CI](https://github.com/light-suzuki/primerblast-oss/actions/workflows/ci.yml/badge.svg)](https://github.com/light-suzuki/primerblast-oss/actions/workflows/ci.yml)
[![Benchmark](https://github.com/light-suzuki/primerblast-oss/actions/workflows/benchmark.yml/badge.svg)](https://github.com/light-suzuki/primerblast-oss/actions/workflows/benchmark.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)

[English](README.md) | **日本語**

植物育種・遺伝学向けの、ローカルかつオープンソースな **Primer-BLAST 相当のコマンドラインツール**。**Primer3** でPCRプライマーを設計し、ローカルのBLAST+データベースに対して**特異性**を完全オフラインで検証します — 未公開ゲノムや複数品種の同時スクリーニングを含みます。さらに in-silico PCR、領域全体のタイリング、プライマー下SNP検出、アンプリコン保存性、CAPS/dCAPSマーカー設計、実験者向けリスク評価も備えます。

> コアは純Python(標準ライブラリのみ)で、`primer3_core` と BLAST+ を外部プロセスとして呼び出します。ユニットテストは**外部ツールもデータも不要**です。

## なぜローカルかつオープンソースか

NCBI Primer-BLAST は優れていますが、**オープンソースではなく**、ホスト型のWebサービスとしてのみ動作します。実務上これは、監査・フォーク・自前ホストができず、手元のデータの隣で動かせず、未公開・エンバーゴ中のゲノムを投入できないことを意味します。

外部サービスへの依存は、パイプラインをその可用性やポリシー(レート制限、メンテナンス、時折の障害)に縛りつけます。ワークフローを**ローカル・オフライン**に保てばその依存が消え、実行が完全に再現可能になります — FASTA・GFF3・VCF・BLAST DB・ツールのバージョンを固定でき、キューもログインも不要です。これは、Webツールがそもそも扱えないケース(ローカル・未公開・複数品種のゲノム)で最も効いてきます。

primerblast-oss は **MITライセンス**なので、誰でも読み、動かし、その上に構築できます。

## なぜ作ったか

NCBI Primer-BLAST は2つのことを行います:(1) Primer3 が候補プライマーを設計し、(2) BLAST が各プライマーをデータベースに照合して**ヒットを予測アンプリコンにペアリング**し、意図しないPCR産物を洗い出します。ローカルの「プライマー + BLAST」スクリプトの多くはプライマー単位のBLASTだけを行い、(2) — 実際にオフターゲット産物を検出する部分 — を欠いています。

`primerblast-oss` は (2) の独立実装です — NCBIそのままのアルゴリズムではなく、BLASTアラインメントに基づくプライミングモデルを用います — そして、ローカルで育種指向のワークフローが必要とする要素に焦点を当てています:

| | NCBI Primer-BLAST | PrimerServer2 | primerblast-oss |
|---|---|---|---|
| Primer3 による設計 | ✅ | ✅ | ✅ |
| BLASTヒットを予測アンプリコンにペアリング | ✅ | ✅ | ✅ |
| いずれのプライマー由来でも F/F, R/R, F/R のオフターゲット産物 | ✅ | ✅ | ✅ |
| 3'末端を考慮したプライミングモデル | ✅ | ✅ (`--use-3-end`) | ✅ |
| **未公開/ローカル**ゲノムでオフライン動作 | 困難 | ✅ | ✅ |
| 熱力学的オフターゲット評価(Tmベース) | ✅ | ✅ (コアモデル) | オプション (primer3-py) |
| プールの**多重**プライマーダイマー *チェック* | — | ✅ | ✅ |
| **多重**互換セット *設計*(ターゲット毎に1ペア) | — | — | ✅ |
| 1回の実行で**複数データベース**をスクリーニング | — | 部分的 | ✅ |
| 貼り付けたプライマーからの **in-silico PCR**(向き非依存) | — | ✅ | ✅ |
| 重複アンプリコンで**領域全体をタイリング** | — | — | ✅ |
| オフターゲットのゲル分離性(サイズ差考慮) | — | — | ✅ |
| 育種アッセイ: GFF3/VCF/CAPS/QTL + リスク | — | — | ✅ |
| スクリプト可能なCLI + ライブラリ、キュー/ログイン不要 | 限定的 | ✅ | ✅ |
| キュレーション済み・継続更新されるデータベース | ✅ | — | — |
| 成熟したホスト型Webサーバー | ✅ | ✅ | ローカルGUIのみ |

これは**ローカル・オフラインのワークフローへの適合**についての話であって、総合的に優れているという主張ではありません。NCBI Primer-BLAST には本ツールにない実利があります:キュレーションされ継続更新されるデータベース、成熟した熱力学モデル、より踏み込んだプライマーダイマー/ヘアピン解析です。[PrimerServer2](https://github.com/billzt/PrimerServer2) も強力なローカルツールでコアのレシピを多く共有しています。primerblast-oss がその上に足しているのは、領域全体のタイリング、ゲル分離性、1回実行のマルチDBスクリーニング、育種アッセイ(GFF3/VCF/CAPS/QTL/risk)、そしてダイマーの*チェック*だけでなく多重互換セットの*設計*です。複数列に ✅ がある場合、その能力が各側に存在することを意味するだけで、内部モデルが同一であるとか出力が一致することを意味しません。

**ベンチマーク(要約):** ランダムに配置した **40座位の Arabidopsis TAIR10** で、primerblast-oss と PrimerServer2 は**非反復座位の92%**で同一のアンプリコン集合(数・サイズ・座標)を予測しました。公開ゲノムの *Lotus japonicus* では手検証3ペアが完全一致し、Arabidopsisの1座位ではトップの新規設計ペアが PrimerServer2 と**ライブのNCBI Primer-BLAST**の両方と一致しました。手法・数値・残差の分析は [`benchmarks/RESULTS.md`](benchmarks/RESULTS.md) §7–§9 にあります。

## 検証状況

短期目標は、ローカルでスクリプト可能なプライマー作業における **PrimerServer2 のスーパーセット**であること:ローカルゲノムでその特異性挙動に一致しつつ、マルチDBスクリーニング、タイリング、マーカー設計、育種アッセイ出力、オフライン再現性を加えること。現時点の根拠:

- **PrimerServer2 との40座位自動対決(Arabidopsis TAIR10)。** パラメータを揃えた条件で、両ツールは**非反復座位 33/36(92%)**で予測アンプリコン集合が完全一致しました。残る不一致はすべて説明可能です:両ツールともプライマーを非特異と判定するが反復コピーの列挙が異なる反復座位と、プライマーの3'末端が完全には整合しないマージナルなサイト(primerblast-oss はこれを非プライミングとして棄却し、PrimerServer2 は二本鎖Tmで保持)です。実装エラーに起因するものはありません([`benchmarks/RESULTS.md`](benchmarks/RESULTS.md) §9)。
- **PrimerServer2、*Lotus japonicus*。** 手検証3ペアがアンプリコン数・サイズ・座標で完全一致(§7)。
- **NCBI Primer-BLAST(単一座位)。** 公開の Arabidopsis TAIR10 テンプレートで、トップの新規設計ペアはライブのNCBIサービスが返すペアと同一で、両者とも特異的と判定しました。これはスポットチェックであって系統的検証ではありません(§8)。
- **継続的な回帰ベンチマーク。** CIが合成のFASTA/BLASTデータベースを構築し、Primer3設計、BLASTアンプリコンペアリング、重複/オフターゲット分類、熱力学ゲート、多重ダイマーチェックを毎プッシュで検査します。

**NCBI Primer-BLAST** に対してはドロップイン等価を主張しません:NCBIはキュレーション済みで継続更新されるデータベース、ホスト型UX、非公開で長年成熟した特異性モデルで優位を保ちます。primerblast-oss は、重要なデータがローカル・未公開・複数リファレンスである、あるいはスクリプトで再現可能に実行する必要がある場合に、より強い選択肢です。

## 特異性の判定方法

各プライマーペア・各データベースについて:

1. `blastn -task blastn-short` が各プライマーのほぼ全長ヒットをすべて検出します。
2. ヒットが**プライミングサイト**になるのは、プライマーの**3'末端がアラインしており**、その**3'末端塩基が一致**し、3'ウィンドウ(`--three-prime-window`、既定5)内のミスマッチが `≤ --max-3prime-mismatch`(既定1)、プライマー全長のミスマッチ合計が `≤ --max-total-mismatch`(既定4)の場合のみです。アラインしていない5'側塩基はミスマッチとして数えます。
3. 各サブジェクト上で、プラス鎖のプライミングサイトは、産物サイズ窓(`--min-product`..`--max-product`)内の下流のマイナス鎖プライミングサイトすべてとペアリングされます。産物長は5'→5'で測定します(真のPCRアンプリコン長)。
4. 産物が**オンターゲット**なのは、両プライマーが完全アニール(ミスマッチ0)、F/Rの向き、設計サイズ(± `--size-tolerance`)のときです。それ以外はすべて**オフターゲット**です。ペアが**特異的**なのは、スクリーニングした全データベースでちょうど1つの産物(意図した産物)だけが予測されるときです。

ペアは、特異性・Tmバランス・GC・3'ダイマー強度でスコア化され A–D にランク付けされます。

`--specificity-profile ncbi` を使うと、ミスマッチ閾値を NCBI Primer-BLAST 風のストリンジェンシープロファイルに切り替えます:合計5ミスマッチまでを候補プライミングサイトとして残し、3'末端5 bp内のミスマッチを1つまで許容し、末端塩基のミスマッチを即棄却せずにカウントします。これは閾値挙動の互換プロファイルであって、NCBIの非公開アルゴリズムやデータベースそのものではありません。

## 必要環境

- Python ≥ 3.8(標準ライブラリのみ)
- `primer3_core`(Debian/Ubuntu: `apt install primer3`)
- BLAST+ の `blastn` / `makeblastdb`(`apt install ncbi-blast+`)
- ヌクレオチドBLASTデータベース(下記参照)

## インストール

```bash
pip install -e .              # `primerblast-oss` コマンドを提供
pip install -e '.[thermo]'    # + オプションの primer3-py(熱力学とダイマー)
# インストールせずに実行する場合:
python -m primerblast_oss --help
```

インストール後の `primerblast-oss <subcommand>` と `python -m primerblast_oss <subcommand>` は等価です。本READMEでは、インストールなしでも例が動くよう `python -m` 形式を使います。

## クイックスタート

```bash
# 1. ゲノムFASTAからBLASTデータベースを構築(初回のみ)
python -m primerblast_oss makedb genome.fa --out-db mydb

# 2. テンプレート上でプライマーを設計し、そのゲノムに対して検証
python -m primerblast_oss design \
  --template-fasta my_gene.fa --product-size 150-500 --db mydb

# 3. または、手持ちのペアを in-silico PCR するだけ
python -m primerblast_oss check \
  --forward GACAAGGAATCAGCGGCTCT --reverse GCAGCGTTTTGTAGTGGGTG --db mydb
```

これらのサブコマンドをラップするローカルのブラウザGUIもありますが、あくまでオプションのおまけです — 末尾の [Web GUI(オプション)](#web-guiオプション) を参照してください。

## 使い方

primerblast-oss は **CLIツール**です。サブコマンド:**design**、**check**、**multiplex**、**multiplex-design**、**tile**、**assay**、**markers**、**makedb**。各サブコマンドの全オプションは `python -m primerblast_oss <subcommand> --help` で確認できます。

`multiplex` はプールされたプライマー間のプライマーダイマー互換性をチェックします(`primer3-py` が必要) — 全プライマー総当たりで、一緒に実行できるセットを選びます:

```bash
python -m primerblast_oss multiplex \
  --primer A_F=... --primer A_R=... --primer B_F=... --primer B_R=...
```

`multiplex-design` はさらに一歩進みます:複数のターゲット(複数レコードのテンプレートFASTA)を渡すと、各ターゲットの候補を設計し、**ターゲット毎に互いに互換な1ペア**を、どの2プライマーも懸念すべき交差ダイマーを作らないように選びます。NCBI Primer-BLAST は各アンプリコンを独立に設計するため、これはできません。

```bash
python -m primerblast_oss multiplex-design \
  --template-fasta targets.fa --db $DB --genome-fasta genome.fa \
  --product-size 80-300 --candidates-per-target 5 --require-specific
```

### `design` — 領域 + 産物サイズ → プライマーペア

```bash
python -m primerblast_oss design \
  --template-fasta my_gene.fa \
  --db /path/to/genome_db \
  --product-size 150-500 --format text
```

複数品種のゲノムを一度にスクリーニング(*すべて*で特異的):

```bash
python -m primerblast_oss design \
  --template "ACGT..." --template-id MyLocus \
  --db /data/blastdb/cultivarA --db /data/blastdb/cultivarB \
  --product-size 200-800 --format tsv
```

### `check` — プライマー配列 → 予測される全PCR産物(in-silico PCR)

プライマーを貼り付けます。向きは**拘束されません**(どのプライマーもフォワード/リバースとして働けます)。各産物はサイズと最近接産物までのサイズ差とともに列挙されるので、余分なバンドが分離可能かを判断できます。

```bash
python -m primerblast_oss check \
  --forward GCACTCTAGAGGTTCAAGGCC --reverse TGGTACGTGTGGTTCAGTTTCA \
  --db /path/to/genome_db
# または名前付きのプライマープール:
python -m primerblast_oss check \
  --primer F1=ACGT... --primer F2=TTGC... --primer R1=GGCA... \
  --db /path/to/genome_db --format json
```

### `tile` — 領域 + アンプリコン長 → 領域全体を覆う重複アンプリコン

Primer-BLAST はターゲット周辺に1つのアンプリコンを設計しますが、`tile` は代わりに領域全体を重複アンプリコンで走査します(例:遺伝子全体のシーケンス)。

```bash
python -m primerblast_oss tile \
  --template-fasta gene.fa \
  --amplicon-min 400 --amplicon-max 700 --overlap 60 \
  --db /path/to/genome_db
```

### `assay` — 遺伝子 / 区間 / SNP からの完全な育種アッセイ

**ローカルゲノム + GFF3/VCF** からターゲットを解決し、プライマーを設計し、**複数のリファレンスゲノム**にわたって特異性をチェックし、**プライマー下のSNP**を洗い出し、**アンプリコン保存性**をスコア化し、オプションの **CAPS/dCAPS** 酵素スキャンを実行し、実験者向けの**リスク(low/medium/high)**を割り当てます。出力は text、JSON、CSV、BED、オリゴの**発注表**、または自己完結型の **HTML** レポートです。

```bash
DB=/path/to/blastdb
# 遺伝子を3品種でスクリーニング、VCFとCDS特徴つき
python -m primerblast_oss assay \
  --gene Psat.cameor.v2.1g00050 --gene-feature cds --gff3 genome.gff3 \
  --genome genome.fa \
  --db $DB/cameor_v2 --db $DB/unpublished_cultivar --db $DB/ZW6 \
  --vcf variants.vcf --flank 100 --product-size 150-600 --format html --out report.html

# SNPをまたぐCAPSマーカー(alt対立遺伝子を指定)
python -m primerblast_oss assay --snp chr1:6385 --alt A \
  --genome genome.fa --db $DB/cameor_v2 --flank 250 --format text
```

### `markers` — QTL区間にわたる等間隔マーカー

```bash
python -m primerblast_oss markers --interval chr1:80000000-90000000 \
  --genome genome.fa --db $DB/cameor_v2 --n-markers 20 --format json
```

### `makedb` — データベース構築(`-parse_seqids` つき)

```bash
python -m primerblast_oss makedb genome.fa --out-db genome_db
```

## NCBI Primer-BLAST のよくある不満点への対処

以下は本ツールが設計上狙っている不満点です。カバレッジには差があり、一部は部分的です — [制限事項](#制限事項) を参照。

| # | よくある不満点 | primerblast-oss のアプローチ |
|---|---|---|
| 1 | バッチ/多数領域に弱い | BED/遺伝子リストにわたる `markers`・`assay`;CLI + ライブラリでスクリプト可能 |
| 2 | ローカル/独自アセンブリへの適合性が低い | すべてローカルFASTA + BLAST DBで動作;`makedb` ヘルパー |
| 3 | 複数リファレンス比較に弱い | `--db` を繰り返し指定可;リファレンス毎にアンプリコン**保存性**をスコア化 |
| 4 | 染色体全体の設計が不便 | `tile` + `markers` が領域/区間全体にプライマーを生成 |
| 5 | プライマーの鎖/向きが不明瞭 | 各結合サイトが鎖・5'/3'座標・伸長方向を報告 |
| 6 | 予期しない副産物が読みにくい | BLASTヒットをサイズ付きの**予測アンプリコンにペアリング** |
| 7 | F-F / R-R 産物が見えにくい | 明示的に列挙しASCIIマップ・表に表示 |
| 8 | 3'末端ミスマッチが見えない | ヒット毎に明示的な3'末端 **5 bp / 10 bp** ミスマッチ数 |
| 9 | パラログ/反復/重複 | ゲノム全体のペアリングが重複プライミングサイトを可視化(専用のリピートマスクなし;BLAST `-max_target_seqs` に律速) |
| 10 | CAPS/dCAPS 向けに作られていない | `caps` スキャン:2対立遺伝子を異なる形で消化する酵素、ゲル差 |
| 11 | GFF3 / VCF / QTL 連携が弱い | `--gene`/`--gff3`、`--vcf`、`--interval`、BED入力 |
| 12 | 空の結果が不透明 | Primer3 の explain 文字列を提示;ステージ毎の診断 |
| 13 | 再現性が弱い | プロベナンスマニフェストがツールバージョン・パラメータ・DBフィンガープリントを固定 |
| 14 | 実験者向けスコアが弱い | `risk` が全シグナルを理由つきで low/medium/high に集約 |
| 15 | 副産物が可視化されない | ASCIIオフターゲットマップ + ゲノムブラウザ用BEDトラック |

設計のヒントは [PrimerServer2](https://github.com/billzt/PrimerServer2)(鎖考慮のBLASTヒットペアリング、マルチスレッド `blastn`、座標入力)と NCBI Primer-BLAST から得ました。

## 制限事項

本ツールが*行わない*ことがわかるよう、正直にスコープを示します:

- **NCBI Primer-BLAST に対してはスポットチェックのみ。** これは独立実装です。公開Arabidopsisの1座位はNCBIと完全一致しました(同一のトップペア、同一の特異性判定 — `benchmarks/RESULTS.md` §8)が、それは単一座位であって系統的な検証ではありません。それ以外の場所の結果はもっともらしく内部チェック済みですが、NCBIの出力と一致する保証はありません。
- **熱力学評価はオプション**(`pip install primer3-py` が必要)。`--genome-fasta` を与えると(`assay` では自動)、各サイトに primer3 による二本鎖Tmと3'末端ΔGが付き、熱力学的に成立しないサイトはゲートで除外されます。無い場合、プライミングはBLASTアラインメントのミスマッチ/3'アンカー規則にフォールバックします。
- **オフターゲット探索は BLAST `-max_target_seqs`(既定5000)に律速。** 極端に反復的な領域では一部ヒットを取りこぼす可能性があり、専用のリピートマスクはありません。
- **dCAPS対応はベストエフォート**で、CAPS判定は酵素テーブル(約40の一般的酵素)に依存し、網羅的なREBASEセットではありません。
- **バッチ/QTLモードは動作しますが大規模ではベンチマークされていません**。各ペアがBLAST検索1回分のコストなので、広いスイープはIO/CPUバウンドです。
- **プライマーダイマー/ヘアピン解析には primer3-py が必要**(オプション)。あれば各ペアにヘアピン/セルフダイマー/交差ダイマーのスコアが付き、`multiplex` サブコマンドがプール全体をチェックします。無い場合はPrimer3の設計時制限のみが適用されます。

これらのいずれへの貢献も歓迎します。

### 主なオプション(design/check/tile 共通)

| オプション | 意味 | 既定値 |
|---|---|---|
| `--product-size` (design) | 1つ以上の範囲、例 `150-500,500-1000` | `70-1000` |
| `--amplicon-min/--amplicon-max/--overlap` (tile) | タイリングの形状 | 400/800/40 |
| `--opt-tm/--min-tm/--max-tm` | プライマーの融解温度窓 | 60/57/63 |
| `--specificity-profile` | ミスマッチプリセット:`local-strict` または `ncbi` | `local-strict` |
| `--max-total-mismatch` | オフターゲットがなおプライミングするのに許容するミスマッチ | 4 |
| `--max-3prime-mismatch` | 3'ウィンドウ内で許容するミスマッチ | 1 |
| `--three-prime-window` | 3'ウィンドウのサイズ | 5 |
| `--max-product` | 増幅可能とみなす最大のオフターゲットアンプリコン | 4000 |
| `--max-target-seqs` | BLASTヒット上限;反復的ゲノムでは上げる | 5000 |
| `--exhaustive` | より高いBLASTヒット上限を使う簡便モード | off |
| `--num-threads` | `blastn` のワーカースレッド数 | 4 |
| `--high-copy-hit-threshold` | リピート感受性の警告を出す生BLAST HSP数 | 10000 |
| `--gel-min-gap` | ゲルで2産物を分離するのに必要なサイズ差(bp) | 50 |
| `--no-3prime-terminal` | 3'末端塩基のミスマッチを許容 | off |
| `--genome-fasta` | design/check/tile で primer3-py の熱力学サイト評価を有効化 | off |
| `--min-anneal-tm` / `--max-3p-dg` | 熱力学的オフターゲットゲートの閾値 | 40 / -5 |
| `--dimer-dg-warn` / `--dimer-tm-warn` | assay/multiplex のプライマーダイマー/ヘアピン警告閾値 | -8 / 45 |
| `--format` | `text` \| `json` \| `tsv`(design のみ) | text |

### 判定とスコアリング

design/tile はペアを **A–D** にランク付けします:**A** = 全データベースで単一産物;**B** = 余分な産物はあるがすべてゲルで分離できるだけサイズが離れている;**C/D** = 意図したバンドと**共泳動**するオフターゲットが1つ以上。共泳動オフターゲットは重くペナルティされ、サイズで分離可能なものは軽微 — 実際の使われ方に合わせています。

## JSON出力(スクリプト向け)

`--format json` は、他ツール(やGUI)へパイプするための安定したスキーマを出力します。各オブジェクトは、消費側が再計算なしに必要とするものを保持します:

- **primer**: `forward`/`reverse`、`tm_f`/`tm_r`、`gc_f`/`gc_r`、`left_start`/`right_start`(0始まりのテンプレート位置)、`product_size`、`penalty`。
- **product**(design/check): `subject`、`start`、`end`、`size`、`orientation`(例 `F/R`、`R/R`)、`fwd_mismatch`+`rev_mismatch`、`nearest_gap`(最近接産物までのbp — ゲル分離性の陰影を駆動)。
- **pair verdict**: `rank`、`score`、`specific_all_db`、`gel_distinguishable`、`total_on_target`/`total_off_target`/`total_comigrating`、および `per_db` 内訳。
- **tile**: `index`、`covers` `[start,end]`、`gap_to_prev`(重複>0 / ギャップ<0)、加えて完全なペアオブジェクト — 領域上にアンプリコントラックを描くのに十分です。

## ライブラリAPI

```python
from primerblast_oss import run_pipeline, in_silico_pcr, design_tiling

# 設計 + 特異性
result = run_pipeline("MyLocus", template_seq, ["/data/db/genome"])
for pair in result.pairs:
    print(pair.forward, pair.reverse, pair.specificity["rank"])

# 任意のプライマーからの in-silico PCR
res = in_silico_pcr({"F": "ACGT...", "R": "TTGC..."}, "/data/db/genome")
for a in res["products"]:
    print(a.size, a.subject, a.orientation)

# 領域全体のタイリング
tiles = design_tiling("gene", template_seq, ["/data/db/genome"],
                      amplicon_min=400, amplicon_max=700, overlap=60)
```

## Web GUI(オプション)

主要なインターフェースはCLIです。利便性のために、各サブコマンドをラップする**ローカルのブラウザフロントエンド**があります — クラウドも第三者Python依存もありません(標準ライブラリの `http.server` の上に構築)。これは主役ではなくおまけです。`primer3_core`、`blastn`、BLASTデータベースがある同じマシン(例:WSL内)で実行してください:

```bash
python -m primerblast_oss.webapp             # http://127.0.0.1:8799 で配信、ブラウザを開く
python -m primerblast_oss.webapp --port 9000 --no-browser
```

design / in-silico PCR / tiling / assay / QTL markers / build DB のタブを提供し、`~/.codex/blast_databases`・`~/blast_databases`・`./databases` 配下のBLASTデータベースを自動検出し、各ジョブをバックグラウンドスレッドで実行し、TSV/CSV/BED/JSONのワンクリックダウンロードを提供します。ヘッダーで English / 日本語 を切り替え可能。ループバック(`127.0.0.1`)にのみバインドし、WSL2では既定のlocalhostフォワーディングでWindowsのブラウザから到達できます。GUIが行うことはすべてCLIからも利用できます。

## テスト

ユニットテストは純Pythonで、**外部ツールもデータも不要**です:

```bash
pip install -e ".[dev]"
pytest                                # またはファイルを直接実行:
python tests/test_specificity.py      # 特異性 / アンプリコンペアリング
python tests/test_integration.py      # 変異、保存性、リスク、CAPS
```

## ベンチマーク

`benchmarks/run_benchmark.py` は、`.fai` インデックス付きゲノムから実領域を抽出し、プライマーを設計し、特異性をスクリーニングします — タイミングとペア毎の予測産物数を報告します。`export PBO_DBDIR=/path/to/blastdb` で自分のデータに向けられます。エンドウゲノムでの実行(設計、複数品種、in-silico PCR、タイリング、完全アッセイ、CAPS)は [`benchmarks/RESULTS.md`](benchmarks/RESULTS.md) にまとめています。

`benchmarks/head_to_head_ps2.py` は [PrimerServer2](https://github.com/billzt/PrimerServer2) との自動一致ベンチマークです:ゲノム上のN個の窓それぞれでペアを設計し、条件を揃えて両ツールの予測アンプリコンを比較します([`benchmarks/RESULTS.md`](benchmarks/RESULTS.md) §9 参照)。

```bash
python benchmarks/head_to_head_ps2.py --genome tair10.fa --db tair10.fa \
  --primertool /path/to/primertool --n-loci 40 --out h2h.json
```

`benchmarks/continuous_benchmark.py` はCIフレンドリーな回帰ベンチマークです:小さな合成FASTA/BLASTデータベースを構築し、Primer3設計、BLASTアンプリコンペアリング、重複/オフターゲット分類、オプションの熱力学ゲート、多重ダイマーチェックを走らせます。

```bash
python benchmarks/continuous_benchmark.py --max-seconds 30
```

## コントリビュート

貢献を歓迎します — [CONTRIBUTING.md](CONTRIBUTING.md) を参照。変更は [CHANGELOG.md](CHANGELOG.md) で追跡しています。

## 引用

研究で利用する場合は引用してください([CITATION.cff](CITATION.cff) 参照)。

## ライセンス

[MIT](LICENSE) © primerblast-oss contributors
