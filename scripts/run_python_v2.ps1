[CmdletBinding()]
param(
    [Parameter(Mandatory = $true, Position = 0)]
    [string]$Script,
    [Alias('Args')]
    [string[]]$Arguments = @(),
    [string]$PythonExe,
    [string]$RepoRoot,
    [ValidateRange(1, 240)]
    [int]$TimeoutMinutes = 30
)

$ErrorActionPreference = 'Stop'

if (-not $RepoRoot) {
    $RepoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
    $RepoRoot = Split-Path -Parent $RepoRoot
}
$RepoRoot = (Resolve-Path -LiteralPath $RepoRoot).Path
$v2Root = Join-Path $RepoRoot 'python_bcu_v2'

if (-not $PythonExe) {
    $venvPython = Join-Path $RepoRoot '.venv\Scripts\python.exe'
    if (Test-Path -LiteralPath $venvPython -PathType Leaf) {
        $PythonExe = $venvPython
    } elseif ($env:BCU_PYTHON_EXE -and
              (Test-Path -LiteralPath $env:BCU_PYTHON_EXE -PathType Leaf)) {
        $PythonExe = $env:BCU_PYTHON_EXE
    } else {
        $candidate = Get-Command python -ErrorAction SilentlyContinue
        if ($candidate) { $PythonExe = $candidate.Source }
    }
}
if (-not $PythonExe -or -not (Test-Path -LiteralPath $PythonExe -PathType Leaf)) {
    throw '找不到 Python。请先创建 .venv，或用 -PythonExe 指定 Python 3.12 的绝对路径。'
}
$PythonExe = (Resolve-Path -LiteralPath $PythonExe).Path

$scriptPath = Join-Path $v2Root $Script
if (-not (Test-Path -LiteralPath $scriptPath -PathType Leaf)) {
    throw "v2 脚本不存在: $scriptPath"
}

# Use ProcessStartInfo instead of Start-Process.  This keeps the explicit
# working directory and avoids the non-ASCII-path startup failure seen on this
# Windows profile.  The v2 scripts also need the sibling python_bcu directory
# discovered from this working directory.
$psi = [System.Diagnostics.ProcessStartInfo]::new()
$psi.FileName = $PythonExe
$psi.WorkingDirectory = $v2Root
$psi.UseShellExecute = $false
$psi.CreateNoWindow = $true
$psi.RedirectStandardOutput = $true
$psi.RedirectStandardError = $true
$psi.EnvironmentVariables['PYTHONUTF8'] = '1'
$psi.EnvironmentVariables['MPLBACKEND'] = 'Agg'

$argsAll = @('-B', $Script) + @($Arguments)
$argListProperty = $psi.PSObject.Properties['ArgumentList']
if ($null -ne $argListProperty) {
    foreach ($arg in $argsAll) { [void]$psi.ArgumentList.Add([string]$arg) }
} else {
    # Windows PowerShell 5.1 has no ArgumentList property.  All current v2
    # entry points use simple scalar arguments; quote each token for paths.
    $psi.Arguments = (($argsAll | ForEach-Object {
        '"' + ([string]$_ -replace '"', '\"') + '"'
    }) -join ' ')
}

$proc = [System.Diagnostics.Process]::new()
$proc.StartInfo = $psi
$started = $false
try {
    $started = $proc.Start()
    if (-not $started) { throw "无法启动 Python: $PythonExe" }
    $stdoutTask = $proc.StandardOutput.ReadToEndAsync()
    $stderrTask = $proc.StandardError.ReadToEndAsync()
    $finished = $proc.WaitForExit([Math]::Max(1, $TimeoutMinutes) * 60 * 1000)
    if (-not $finished) {
        try { $proc.Kill() } catch { }
        Write-Error "Python 运行超过 $TimeoutMinutes 分钟，状态应记录为 BLOCKED。"
        exit 124
    }
    $stdoutTask.Wait(); $stderrTask.Wait()
    if ($stdoutTask.Result) { [Console]::Out.Write($stdoutTask.Result) }
    if ($stderrTask.Result) { [Console]::Error.Write($stderrTask.Result) }
    exit $proc.ExitCode
} finally {
    $proc.Dispose()
}
