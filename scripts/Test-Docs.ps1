[CmdletBinding()]
param(
    [string]$RepoRoot,
    [int]$ExpectedSourceFiles = 143
)

$ErrorActionPreference = 'Stop'
if (-not $RepoRoot) {
    $RepoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
    $RepoRoot = Split-Path -Parent $RepoRoot
}
$failures = [System.Collections.Generic.List[string]]::new()
$root = (Resolve-Path -LiteralPath $RepoRoot).Path
function Fail([string]$Message) { $failures.Add($Message) }

$required = @(
    'README.md', 'SOURCE_MANIFEST.csv',
    'docs/01_项目现状与边界.md', 'docs/02_开发交接.md', 'docs/03_验证覆盖矩阵.md',
    'docs/04_操作与回归指南.md', 'docs/05_开发路线图.md', 'docs/06_实验与结果规范.md',
    'docs/修改日志.md', 'docs/provenance/2026-09-01_接手来源清单.md',
    'scripts/Import-SourceDocs.ps1', 'scripts/audit_environment.ps1',
    'scripts/bootstrap_local.ps1', 'scripts/run_all_validation.ps1'
)
foreach ($relative in $required) {
    if (-not (Test-Path -LiteralPath (Join-Path $root $relative) -PathType Leaf)) { Fail "missing required file: $relative" }
}

$manifestPath = Join-Path $root 'SOURCE_MANIFEST.csv'
$rows = @()
if (Test-Path -LiteralPath $manifestPath) { $rows = @(Import-Csv -LiteralPath $manifestPath) }
if ($rows.Count -ne $ExpectedSourceFiles) { Fail "manifest row count $($rows.Count) != expected $ExpectedSourceFiles" }

function Get-LineCount([string]$Path) {
    $bytes = [System.IO.File]::ReadAllBytes($Path)
    if ($bytes.Length -eq 0) { return 0 }
    $n = 0; foreach ($b in $bytes) { if ($b -eq 10) { $n++ } }
    if ($bytes[$bytes.Length - 1] -ne 10) { $n++ }
    return $n
}
foreach ($row in $rows) {
    $destination = Join-Path $root ($row.destination_path -replace '/', '\')
    if (-not (Test-Path -LiteralPath $destination -PathType Leaf)) { Fail "manifest destination missing: $($row.destination_path)"; continue }
    $item = Get-Item -LiteralPath $destination
    $hash = (Get-FileHash -LiteralPath $destination -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($hash -ne $row.sha256.ToLowerInvariant()) {
        # python_bcu_v2 is the explicitly designated development platform in
        # this integrated repository.  Its baseline hash remains in the
        # manifest for provenance, while intentional working-tree changes are
        # reviewed by Git.  MATLAB and v1 snapshots stay immutable here.
        if ($row.destination_path -notlike 'python_bcu_v2/*') {
            Fail "SHA-256 mismatch: $($row.destination_path)"
        }
    }
    if ($row.destination_path -notlike 'python_bcu_v2/*') {
        if ([int64]$row.bytes -ne $item.Length) { Fail "byte count mismatch: $($row.destination_path)" }
        if ([int]$row.lines -ne (Get-LineCount $destination)) { Fail "line count mismatch: $($row.destination_path)" }
    }
}

$allowedExtensions = @('.md','.py','.m','.ps1','.json','.csv','.yaml','.yml','.toml','.txt','.ps','.gitignore','')
$allowedTopLevel = @(
    'README.md', 'SOURCE_MANIFEST.csv', '.gitignore',
    'docs', 'sources', 'scripts', 'matlab_platform', 'python_bcu',
    'python_bcu_v2', 'validation', 'experiments'
)
# The integrated checkout may be created by the sandbox service account while
# validation runs as the user's account.  Pass an explicit, narrow safe
# directory so a dubious-ownership warning cannot silently turn the tracked
# file list into an empty list.
$tracked = @(git -c "safe.directory=$root" -c core.quotePath=false -C $root ls-files)
foreach ($relative in $tracked) {
    $normalized = $relative -replace '\\', '/'
    $parts = $normalized.Split('/', 2)
    if ($parts.Count -eq 1) {
        if ($allowedTopLevel -notcontains $parts[0]) {
            Fail "unapproved tracked root file: $relative"
        }
    } elseif ($allowedTopLevel -notcontains $parts[0]) {
        Fail "unapproved tracked top-level directory: $relative"
    }
    if ($normalized -match '(^|/)(C1_Matpower|results|figures|__pycache__|\.git)(/|$)' -or
        $normalized -match '\.(mat|pyc|png|jpg|jpeg|gif|pdf|zip)$') {
        Fail "forbidden tracked path: $relative"
    }
    $extension = [System.IO.Path]::GetExtension($relative).ToLowerInvariant()
    if ($extension -eq '' -and [System.IO.Path]::GetFileName($relative) -ne '.gitignore') { $extension = '' }
    if ($allowedExtensions -notcontains $extension -and [System.IO.Path]::GetFileName($relative) -ne '.gitignore') {
        Fail "unapproved tracked extension: $relative"
    }
}

$outdated = @(
    '目前只有 `reduced_cct` 一条路径验证过',
    'SPM 能量法在 Python 里根本没实现',
    'reduced region、reduced numerical、two-machine 尚未对照',
    'v1 仍绝对冻结'
)
$authoritative = @('README.md')
$authoritative += @(Get-ChildItem -LiteralPath (Join-Path $root 'docs') -Recurse -File -Filter '*.md' |
    ForEach-Object { $_.FullName.Substring($root.Length + 1) })
foreach ($relative in $authoritative) {
    $path = Join-Path $root $relative
    if (-not (Test-Path -LiteralPath $path)) { continue }
    $text = Get-Content -LiteralPath $path -Raw -Encoding utf8
    foreach ($phrase in $outdated) { if ($text.Contains($phrase)) { Fail "outdated authoritative wording in $relative : $phrase" } }
}

# 检查仓库内部 Markdown 链接；外部 URL、锚点和代码片段不参与。
$markdownFiles = Get-ChildItem -LiteralPath $root -Recurse -File -Filter '*.md' |
    Where-Object {
        $_.FullName -notmatch '[\\/]C1_Matpower([\\/]|$)' -and
        $_.FullName -notmatch '[\\/](results|figures|__pycache__)([\\/]|$)'
    }
foreach ($file in $markdownFiles) {
    $text = Get-Content -LiteralPath $file.FullName -Raw -Encoding utf8
    foreach ($match in [regex]::Matches($text, '\[[^\]]*\]\(([^)]+)\)')) {
        $target = $match.Groups[1].Value.Trim().Trim('<','>')
        if ($target -match '^(?:https?:|mailto:|#)') { continue }
        if ($target -match '\s' -and $target -notmatch '[\\/]') { continue }
        $targetPath = ($target -split '#', 2)[0]
        if (-not $targetPath) { continue }
        $resolved = [System.IO.Path]::GetFullPath((Join-Path $file.DirectoryName $targetPath))
        if (-not (Test-Path -LiteralPath $resolved -PathType Leaf)) {
            # 来源平台文档有少量相对统一工作树根目录的链接（例如 python_bcu/...）。
            $resolved = [System.IO.Path]::GetFullPath((Join-Path $root $targetPath))
        }
        if (-not (Test-Path -LiteralPath $resolved -PathType Leaf)) {
            Fail "broken Markdown link in $($file.FullName.Substring($root.Length + 1)): $target"
        }
    }
}

if ($failures.Count -gt 0) {
    $failures | ForEach-Object { Write-Error $_ }
    exit 1
}
Write-Output "DOCS_OK manifest=$($rows.Count) tracked=$($tracked.Count)"
