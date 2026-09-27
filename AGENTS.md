# AGENTS.md

このファイルは、Codex/自動エージェントが xTB_Analytical_Hessian リポジトリで作業するときの
Python 実行環境、標準コマンド、サブエージェント運用を共有するためのメモです。

## 1. 実行環境（前提）

- OS / shell: Windows + PowerShell
- リポジトリルート: `C:\GitHub_repository\xTB_Analytical_Hessian`
- Codex Agent environment: Windows native
- Python環境: Mambaforge / Conda `xtb` (`C:\Users\f3r1i\mambaforge\envs\xtb\python.exe`)
- Python標準バージョン: 3.13.15

## 2. Pythonコマンドの扱い

- Python は Windows版 Mambaforge / Conda の仮想環境 `xtb` を使用する。
- Python バージョンは `3.13.15` を標準とする。
- 実行ファイルは `C:\Users\f3r1i\mambaforge\envs\xtb\python.exe` を使う。
- Codex / 自動エージェントは `conda activate xtb` の成否や現在の shell activation 状態に依存してはならない。
- Python を実行するときは、可能な限り上記 `python.exe` を直接指定する。
- PowerShell では次の変数を標準として使う。
  - `$PY = "C:\Users\f3r1i\mambaforge\envs\xtb\python.exe"`
- pip は単独の `pip` コマンドではなく、必ず `& $PY -m pip ...` の形で実行する。
- pytest などの Python モジュールも、必ず `& $PY -m <module>` の形で実行する。
- 新しいスクリプトで子プロセス Python を起動する場合は、`sys.executable` を優先する。
- 使用中の Python を確認するときは `& $PY -c "import sys; print(sys.executable)"` を使い、上記パスと一致することを確認する。
- system Python、Microsoft Store Python、WSL Python、WinPython、他の Conda 環境を xTB_Analytical_Hessian の開発・テストに使用しない。

## 2.5 環境確認

作業開始時、Python実行や依存関係の変更を伴う場合は、必要に応じて次を確認する。

```powershell
$PY = "C:\Users\f3r1i\mambaforge\envs\xtb\python.exe"
& $PY --version
& $PY -c "import sys; print(sys.executable)"
```

期待する実行ファイル:

```text
C:\Users\f3r1i\mambaforge\envs\xtb\python.exe
```

## 3. 依存関係

- Python パッケージメタデータと通常の開発依存関係の正本は `pyproject.toml` とする。
- 開発環境への標準インストールは、リポジトリルートで次を使う。
  - `& $PY -m pip install -e ".[dev]"`
- 既存の requirements / lock ファイルを利用する処理がある場合は、その用途を維持する。
- Windows配布や厳密な再現性が必要な場合は、リポジトリ内の lock ファイルを優先する。
- 依存関係を追加・変更する場合は、まず `pyproject.toml` を更新し、必要に応じて lock ファイルも更新する。
- Codex / 自動エージェントは、依存関係の導入先が `xtb` 環境であることを確認してから install / update を行う。

## 3.5 正本ドキュメント

実装前に `docs/index.md` を入口として必要な範囲を参照する。AI専用のarchitecture/test仕様を別ファイルへ複製しない。

- ドキュメント入口: `docs/index.md`
- 自動エージェント入口: `docs/developer/codex.md`
- 製品仕様: `docs/reference/specification.md`
- architecture正本: `docs/developer/architecture.md`
- テスト/チェックのルーティング: `docs/developer/testing.md`
- 開発方針: `docs/developer/development.md`
- リリース運用: `docs/developer/release/release_checklist.md`

`docs/archive/` は過去の計画・計測・受け入れ記録であり、現行仕様の正本として扱わない。

## 4. 標準コマンド

### アプリ起動

```powershell
$PY = "C:\Users\f3r1i\mambaforge\envs\xtb\python.exe"
& $PY main.py
```

### 変更時の基本検証

```powershell
& $PY -m compileall -q .
& $PY -m pytest -q
```

## 4.5 変更時の必須運用ルール

- 単一エージェント作業では、コード、ドキュメント、テスト、パッケージング設定など、リポジトリに変更を加える論理タスクの完了時に `& $PY tools/bump_version.py` を使い、現在のrelease channelを維持したまま `MAJOR.MINOR.PATCH` のpatch番号を 1 つ上げる。
- サブエージェントを使用する場合、planner / implementer / tester は `tools/bump_version.py` を実行しない。複数エージェントの変更を統合した後、親agentだけが論理タスクにつき1回実行する。
- 変更後は、変更範囲に応じた最小検証を先に行う。全標準検証は、論理的な一連の変更が統合された後に tester または親agentが原則1回だけ実行する。
- 同じコミット候補に対して、複数のサブエージェントが同じフルテストを重複実行しない。

## 4.6 サブエージェント運用（省リソース）

Codexは、タスクの複雑さに応じて必要な場合だけサブエージェントを使用する。
小さく明確な変更では、サブエージェントを無理に生成せず、親agentまたは単一workerで完了してよい。

標準ロール:

- `planner`: 実装前の調査と必要最小限の計画。コードは変更しない。
- `implementer`: 合意済みの範囲だけを実装する。
- `tester`: テスト実行、失敗の切り分け、簡潔な結果報告を行う。原則としてコードを変更しない。

通常フロー:

```text
小変更:
parent -> implementer -> 必要な局所テスト

中規模以上:
parent -> planner -> implementer -> tester -> parent final verification

難しい不具合:
上記フローで解決できない場合だけ、より高性能なモデル/高いreasoningへ昇格
```

リソース節約ルール:

- 同じ調査を複数agentに重複させない。
- plannerが特定した変更対象をimplementerがゼロから再探索しない。
- implementerは変更箇所に直接関係する最小テストだけを実行する。
- testerは最初に狭いテストを実行し、失敗や影響範囲が必要な場合だけテスト範囲を広げる。
- テスト成功時は長いstdout/stderrを親agentへ返さない。実行コマンド、exit code、pass/fail/skipの要約だけ返す。
- テスト失敗時は、failing test名、重要なエラー抜粋、推定原因だけを返す。
- リポジトリ全体のフルテストは、同一変更セットについて原則1回だけ実行する。
- 不要なレビュー専用agent、重複したexplorer、過剰な並列化を生成しない。
- サブエージェント同士でさらにサブエージェントを増殖させない。必要な追加委譲は親agentへ返す。
- 明確な作業は低コストモデル/低reasoningを優先し、難易度が上がった場合だけ段階的に昇格する。
- 高性能モデルへの昇格条件は、architecture横断変更、再現困難な不具合、複数回の合理的な修正でも解決しない問題、または重要な設計判断に限定する。

モデルの標準割り当ては `.codex/agents/*.toml` を正本とする。
親agentのモデルはCodex App/セッション側で選択し、プロジェクト設定では固定しない。


## 4.7 デフォルト実装ワークフロー

コード変更を伴う依頼では、ユーザーから別の進め方が明示されていない限り、以下をデフォルトとする。

- 小さく明確な変更:
  - `planner` は省略してよい。
  - `implementer` が実装し、変更箇所に直接関係する必要最小限の局所テストを行う。
  - 変更範囲が十分に限定されており、追加の独立検証が不要な場合は `tester` も省略してよい。
- 中規模以上の変更:
  - 原則として `planner -> implementer -> tester` の順で進める。
  - `planner` が変更対象・実装手順・必要テストを特定した後、`implementer` はその計画を再調査せず実装へ進む。
  - `tester` は統合後の変更セットに対し、狭いテストから開始して必要な場合だけ範囲を広げる。
- `planner`、`implementer`、`tester` は `.codex/agents/*.toml` に定義されたモデル、reasoning effort、sandbox設定を使用する。
- 同じ調査、実装、テストを複数agentで重複させない。
- 親agentは、サブエージェントで十分に完了した作業をゼロからやり直さない。
- サブエージェントは追加のサブエージェントを生成しない。追加調査や高性能モデルが必要な場合は親agentへ返す。
- Terra で十分に解決できない場合のみ親agentへエスカレーションし、より高性能なモデルまたは高い reasoning effort への昇格を検討する。
- 高性能モデルへの昇格は、architecture横断変更、再現困難な不具合、複数回の合理的な修正でも解決しない問題、重要な設計判断に限定する。
- version bump はサブエージェントでは実行せず、論理タスクの完了時に親agentが原則1回だけ実行する。
- フル標準検証は同一変更セットに対して原則1回だけ実行する。
- ユーザーが「計画だけ」「実装だけ」「テストだけ」など工程を限定した場合は、その指定を優先し、不要なagentを起動しない。

