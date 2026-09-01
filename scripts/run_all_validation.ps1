[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$PythonExe,
    [int]$TimeoutMinutes = 30,
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot)
)

$ErrorActionPreference = 'Stop'
if (-not (Test-Path -LiteralPath $PythonExe -PathType Leaf)) { throw "Python executable does not exist: $PythonExe" }
$stamp = (Get-Date).ToUniversalTime().ToString('yyyyMMdd_HHmmssZ')
$reportDir = Join-Path $RepoRoot "validation/reports/$stamp"
New-Item -ItemType Directory -Path $reportDir -Force | Out-Null
$tasks = @(
    @{ Path = 'python_bcu/tests/smoke_test.py'; Cwd = 'python_bcu' },
    @{ Path = 'python_bcu_v2/run_validation.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/run_matlab_xval.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/test_fixes.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/test_spm_dae.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/test_solvers.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/test_p2.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/test_config.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/test_cuep.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/run_full_xval.py'; Cwd = 'python_bcu_v2' }
)
$records = [System.Collections.Generic.List[object]]::new()
foreach ($task in $tasks) {
    $scriptPath = Join-Path $RepoRoot $task.Path
    $stdout = Join-Path $reportDir ((Split-Path $task.Path -Leaf) + '.stdout.txt')
    $stderr = Join-Path $reportDir ((Split-Path $task.Path -Leaf) + '.stderr.txt')
    $start = Get-Date
    $status = 'FAILED'; $exitCode = -1
    if (-not (Test-Path -LiteralPath $scriptPath -PathType Leaf)) {
        Set-Content -LiteralPath $stderr -Value 'script not found' -Encoding utf8
    } else {
        $proc = Start-Process -FilePath $PythonExe -ArgumentList @('-B', (Split-Path $task.Path -Leaf)) -WorkingDirectory (Join-Path $RepoRoot $task.Cwd) -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
        $finished = $proc.WaitForExit([Math]::Max(1, $TimeoutMinutes) * 60 * 1000)
        if (-not $finished) {
            $proc.Kill()
            $status = 'BLOCKED'
            $exitCode = -2
        } else {
            $exitCode = $proc.ExitCode
            $status = if ($exitCode -eq 0) { 'PASSED' } else { 'FAILED' }
        }
    }
    $end = Get-Date
    $records.Add([pscustomobject]@{
        command = "$PythonExe -B $($task.Path)"
        working_directory = (Join-Path $RepoRoot $task.Cwd)
        start_time = $start.ToUniversalTime().ToString('o')
        end_time = $end.ToUniversalTime().ToString('o')
        exit_code = $exitCode
        status = $status
        stdout_path = $stdout.Substring($RepoRoot.Length + 1)
        stderr_path = $stderr.Substring($RepoRoot.Length + 1)
    })
}
$records | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $reportDir 'summary.json') -Encoding utf8
$records | Format-Table -AutoSize
if (@($records | Where-Object { $_.status -in @('FAILED','BLOCKED') }).Count -gt 0) { exit 1 }
