[CmdletBinding()]
param(
    [string]$PythonExe,
    [string]$MatlabExe,
    [string]$MatpowerRoot,
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot)
)

$ErrorActionPreference = 'Stop'
$blocked = $false
$report = [ordered]@{
    audited_at = (Get-Date).ToUniversalTime().ToString('o')
    repo_root = (Resolve-Path -LiteralPath $RepoRoot).Path
    repo_head = $null
    repo_branch = $null
    repo_status = @()
    python = [ordered]@{ executable = $PythonExe; available = $false; version = $null; packages = @{} }
    matlab = [ordered]@{ executable = $MatlabExe; available = $false; version = $null }
    matpower_root = $MatpowerRoot
    local_mat_files = @()
    source_commits = @{}
    status = 'OK'
    blockers = @()
}

if (Test-Path -LiteralPath (Join-Path $RepoRoot '.git')) {
    $report.repo_head = (git -C $RepoRoot rev-parse HEAD 2>$null)
    $report.repo_branch = (git -C $RepoRoot branch --show-current 2>$null)
    $report.repo_status = @(git -C $RepoRoot status --short 2>$null)
}

if (-not $PythonExe) {
    $candidate = Get-Command python -ErrorAction SilentlyContinue
    if ($candidate) { $PythonExe = $candidate.Source }
}
if ($PythonExe -and (Test-Path -LiteralPath $PythonExe -PathType Leaf)) {
    $report.python.executable = (Resolve-Path -LiteralPath $PythonExe).Path
    try {
        $version = (& $PythonExe --version 2>&1 | Out-String).Trim()
        if ($LASTEXITCODE -eq 0) {
            $report.python.available = $true
            $report.python.version = $version
            foreach ($pkg in @('numpy','scipy','matplotlib')) {
                $out = (& $PythonExe -c "import importlib.metadata as m; print(m.version('$pkg'))" 2>&1 | Out-String).Trim()
                $report.python.packages[$pkg] = if ($LASTEXITCODE -eq 0) { $out } else { $null }
            }
        }
    } catch { $report.python.version = $_.Exception.Message }
}
if (-not $report.python.available) {
    $blocked = $true
    $report.blockers += 'Python 无法启动或路径不存在'
}

if (-not $MatlabExe) {
    $candidate = Get-Command matlab -ErrorAction SilentlyContinue
    if ($candidate) { $MatlabExe = $candidate.Source }
}
if ($MatlabExe -and (Test-Path -LiteralPath $MatlabExe -PathType Leaf)) {
    $report.matlab.executable = (Resolve-Path -LiteralPath $MatlabExe).Path
    try {
        $version = (& $MatlabExe -batch "disp(version)" 2>&1 | Out-String).Trim()
        if ($LASTEXITCODE -eq 0) {
            $report.matlab.available = $true
            $report.matlab.version = $version
        }
    } catch { $report.matlab.version = $_.Exception.Message }
}
if (-not $report.matlab.available) {
    $blocked = $true
    $report.blockers += 'MATLAB 无法启动或路径不存在'
}

if ($MatpowerRoot) {
    $report.matpower_root = (Resolve-Path -LiteralPath $MatpowerRoot -ErrorAction SilentlyContinue).Path
    if (-not $report.matpower_root) {
        $blocked = $true
        $report.blockers += 'MATPOWER 路径不存在'
    }
}
$report.local_mat_files = @(Get-ChildItem -LiteralPath $RepoRoot -Recurse -File -Filter '*.mat' -ErrorAction SilentlyContinue |
    ForEach-Object { $_.FullName.Substring($report.repo_root.Length + 1) })
foreach ($d in @('matlab_platform','python_bcu','python_bcu_v2')) {
    $p = Join-Path $RepoRoot $d
    if (Test-Path -LiteralPath (Join-Path $p '.git')) {
        $report.source_commits[$d] = (git -c "safe.directory=$p" -C $p rev-parse HEAD 2>$null)
    }
}
if ($report.source_commits.Count -eq 0) {
    $manifest = Join-Path $RepoRoot 'SOURCE_MANIFEST.csv'
    if (Test-Path -LiteralPath $manifest) {
        foreach ($row in (Import-Csv -LiteralPath $manifest | Select-Object -First 1000)) {
            if (-not $report.source_commits.Contains($row.platform)) {
                $report.source_commits[$row.platform] = $row.source_commit
            }
        }
    }
}

if ($blocked) { $report.status = 'BLOCKED' }
$report | ConvertTo-Json -Depth 8
if ($blocked) { exit 2 }
