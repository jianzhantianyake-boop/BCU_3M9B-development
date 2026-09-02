[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$SourceRoot,

    [string]$RepoRoot,

    [string]$SnapshotDate = (Get-Date -Format 'yyyy-MM-dd')
)

$ErrorActionPreference = 'Stop'
$OutputEncoding = [System.Text.UTF8Encoding]::new()
if (-not $RepoRoot) {
    $RepoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
    $RepoRoot = Split-Path -Parent $RepoRoot
}

if (-not (Test-Path -LiteralPath $SourceRoot -PathType Container)) {
    throw "SourceRoot does not exist: $SourceRoot"
}
if (-not (Test-Path -LiteralPath $RepoRoot -PathType Container)) {
    throw "RepoRoot does not exist: $RepoRoot"
}

$platforms = @(
    [pscustomobject]@{ Name = 'matlab_platform'; Platform = 'matlab_platform' },
    [pscustomobject]@{ Name = 'python_bcu'; Platform = 'python_bcu' },
    [pscustomobject]@{ Name = 'python_bcu_v2'; Platform = 'python_bcu_v2' }
)
$rows = [System.Collections.Generic.List[object]]::new()
$importedAt = (Get-Date).ToUniversalTime().ToString('o')

function Get-LineCount([string]$Path) {
    $bytes = [System.IO.File]::ReadAllBytes($Path)
    if ($bytes.Length -eq 0) { return 0 }
    $count = 0
    foreach ($b in $bytes) { if ($b -eq 10) { $count++ } }
    if ($bytes[$bytes.Length - 1] -ne 10) { $count++ }
    return $count
}

foreach ($platform in $platforms) {
    $sourceRepo = Join-Path $SourceRoot $platform.Name
    $destinationRepo = Join-Path $RepoRoot $platform.Name
    if (-not (Test-Path -LiteralPath $sourceRepo -PathType Container)) {
        throw "Missing source platform: $sourceRepo"
    }

    $sourceCommit = (git -c "safe.directory=$sourceRepo" -C $sourceRepo rev-parse HEAD).Trim()
    if ($LASTEXITCODE -ne 0) { throw "Cannot resolve HEAD for $sourceRepo" }
    $statusLines = @(git -c "safe.directory=$sourceRepo" -C $sourceRepo status --short)
    if ($LASTEXITCODE -ne 0) { throw "Cannot inspect worktree for $sourceRepo" }
    $worktreeState = if ($statusLines.Count -eq 0) { 'clean' } else {
        ('dirty:' + (($statusLines -join '; ') -replace '[\r\n]+', ' '))
    }
    $tracked = @(git -c "safe.directory=$sourceRepo" -c core.quotePath=false -C $sourceRepo ls-files)
    if ($LASTEXITCODE -ne 0) { throw "Cannot enumerate tracked files for $sourceRepo" }

    foreach ($relative in $tracked) {
        if ([string]::IsNullOrWhiteSpace($relative)) { continue }
        $relative = $relative -replace '/', '\'
        if ($relative -match '^(?:C1_Matpower|\.git|results|figures)(?:\\|$)' -or
            $relative -match '\.mat$') { continue }

        $sourcePath = Join-Path $sourceRepo $relative
        $destinationPath = Join-Path $destinationRepo $relative
        if (-not (Test-Path -LiteralPath $sourcePath -PathType Leaf)) {
            throw "Tracked source file is missing: $sourcePath"
        }
        $destinationParent = Split-Path -Parent $destinationPath
        if (-not (Test-Path -LiteralPath $destinationParent)) {
            New-Item -ItemType Directory -Path $destinationParent -Force | Out-Null
        }
        if (Test-Path -LiteralPath $destinationPath -PathType Leaf) {
            $sourceHash = (Get-FileHash -LiteralPath $sourcePath -Algorithm SHA256).Hash
            $destinationHash = (Get-FileHash -LiteralPath $destinationPath -Algorithm SHA256).Hash
            if ($sourceHash -ne $destinationHash) {
                throw "Refusing to overwrite changed destination: $destinationPath"
            }
        } else {
            Copy-Item -LiteralPath $sourcePath -Destination $destinationPath
        }

        $item = Get-Item -LiteralPath $destinationPath
        $repoRelativeSource = ($platform.Name + '\' + $relative) -replace '\\', '/'
        $repoRelativeDestination = $repoRelativeSource
        $rows.Add([pscustomobject]@{
            source_path = $repoRelativeSource
            destination_path = $repoRelativeDestination
            platform = $platform.Platform
            source_commit = $sourceCommit
            source_worktree_state = $worktreeState
            bytes = [int64]$item.Length
            lines = [int](Get-LineCount $destinationPath)
            sha256 = (Get-FileHash -LiteralPath $destinationPath -Algorithm SHA256).Hash.ToLowerInvariant()
            imported_at = $importedAt
        })
    }
}

$manifestPath = Join-Path $RepoRoot 'SOURCE_MANIFEST.csv'
$rows | Sort-Object source_path | ConvertTo-Csv -NoTypeInformation | Set-Content -LiteralPath $manifestPath -Encoding utf8
Write-Output "Imported $($rows.Count) files; manifest: $manifestPath"
