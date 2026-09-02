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
    @{ Path = 'python_bcu/tests/smoke_test.py'; Arg = 'tests/smoke_test.py'; Cwd = 'python_bcu' },
    @{ Path = 'python_bcu_v2/run_validation.py'; Arg = 'run_validation.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/run_matlab_xval.py'; Arg = 'run_matlab_xval.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/test_fixes.py'; Arg = 'test_fixes.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/test_spm_dae.py'; Arg = 'test_spm_dae.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/test_solvers.py'; Arg = 'test_solvers.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/test_p2.py'; Arg = 'test_p2.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/test_config.py'; Arg = 'test_config.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/test_cuep.py'; Arg = 'test_cuep.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/run_full_xval.py'; Arg = 'run_full_xval.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/test_mutations.py'; Arg = 'test_mutations.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'python_bcu_v2/test_spm_xval_gate.py'; Arg = 'test_spm_xval_gate.py'; Cwd = 'python_bcu_v2' },
    @{ Path = 'experiments/test_damping_study.py'; Arg = 'test_damping_study.py'; Cwd = 'experiments' }
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
        # Use ProcessStartInfo with shell execution disabled.  Start-Process
        # can return 0xC0000142 for an otherwise runnable Python executable
        # when this repository is under a non-ASCII Windows path; direct
        # invocation and this API use the same interpreter and working dir.
        $psi = [System.Diagnostics.ProcessStartInfo]::new()
        $psi.FileName = $PythonExe
        $psi.Arguments = '-B "' + $task.Arg + '"'
        $psi.WorkingDirectory = (Join-Path $RepoRoot $task.Cwd)
        $psi.UseShellExecute = $false
        $psi.CreateNoWindow = $true
        $psi.RedirectStandardOutput = $true
        $psi.RedirectStandardError = $true
        $proc = [System.Diagnostics.Process]::new()
        $proc.StartInfo = $psi
        try {
            [void]$proc.Start()
            $stdoutTask = $proc.StandardOutput.ReadToEndAsync()
            $stderrTask = $proc.StandardError.ReadToEndAsync()
            $finished = $proc.WaitForExit([Math]::Max(1, $TimeoutMinutes) * 60 * 1000)
            if (-not $finished) {
                $proc.Kill()
                $status = 'BLOCKED'
                $exitCode = -2
            } else {
                # Wait for redirected streams after process termination so no
                # diagnostic output is lost or deadlocked in a pipe.
                $stdoutTask.Wait(); $stderrTask.Wait()
                Set-Content -LiteralPath $stdout -Value $stdoutTask.Result -Encoding utf8
                Set-Content -LiteralPath $stderr -Value $stderrTask.Result -Encoding utf8
                $exitCode = $proc.ExitCode
                $status = if ($exitCode -eq 0) { 'PASSED' } else { 'FAILED' }
            }
        } catch {
            Set-Content -LiteralPath $stderr -Value $_.Exception.ToString() -Encoding utf8
            $status = 'BLOCKED'
            $exitCode = -2
        } finally {
            $proc.Dispose()
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
