[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$MatpowerRoot,
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot)
)

$ErrorActionPreference = 'Stop'
if (-not (Test-Path -LiteralPath $MatpowerRoot -PathType Container)) {
    throw "MATPOWER source path does not exist: $MatpowerRoot"
}
$destination = Join-Path $RepoRoot 'matlab_platform/C1_Matpower/matpower7.1'
if (Test-Path -LiteralPath $destination) {
    $existing = (Resolve-Path -LiteralPath $destination).Path
    $source = (Resolve-Path -LiteralPath $MatpowerRoot).Path
    if ($existing -ne $source) {
        throw "Local MATPOWER destination already exists and is not the requested source: $existing"
    }
} else {
    New-Item -ItemType Directory -Path (Split-Path -Parent $destination) -Force | Out-Null
    Copy-Item -LiteralPath $MatpowerRoot -Destination $destination -Recurse
}

$required = @('lib','data','most/lib','mp-opt-model/lib','mips/lib','mptest/lib')
$missing = @($required | Where-Object { -not (Test-Path -LiteralPath (Join-Path $destination $_) -PathType Container) })
if ($missing.Count -gt 0) { throw "MATPOWER installation is incomplete: $($missing -join ', ')" }
Write-Output "MATPOWER ready: $destination"
