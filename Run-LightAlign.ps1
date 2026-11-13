[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)][string]$InputFile,
    [Parameter(Mandatory = $true)][string]$OutputDirectory,
    [string]$Executable = '',
    [ValidateRange(1, 4096)][int]$Threads = 12,
    [ValidateRange(1, 2147483647)][int]$GroupSize = 10200,
    [ValidateRange(1, 2147483647)][int]$Window = 30,
    [ValidateRange(0.000001, 0.999999)][double]$ErrorRate = 0.02,
    [ValidateRange(1, 2147483647)][int]$DBA = 87,
    [ValidateRange(1, 2147483647)][int]$UnitLength = 840,
    [ValidateRange(100, 60000)][int]$SampleIntervalMilliseconds = 500
)

Set-StrictMode -Version 2.0
$ErrorActionPreference = 'Stop'
$invariant = [Globalization.CultureInfo]::InvariantCulture
$utf8 = New-Object Text.UTF8Encoding($false)
if ([string]::IsNullOrWhiteSpace($Executable)) {
    $Executable = Join-Path $PSScriptRoot 'bin\LightAlign.exe'
}

function ConvertTo-NativeArgument {
    param([AllowEmptyString()][string]$Value)
    $builder = New-Object Text.StringBuilder
    [void]$builder.Append([char]34)
    $slashes = 0
    foreach ($character in $Value.ToCharArray()) {
        if ($character -eq [char]92) {
            $slashes++
        } elseif ($character -eq [char]34) {
            [void]$builder.Append([char]92, (2 * $slashes + 1))
            [void]$builder.Append([char]34)
            $slashes = 0
        } else {
            if ($slashes -gt 0) { [void]$builder.Append([char]92, $slashes) }
            [void]$builder.Append($character)
            $slashes = 0
        }
    }
    if ($slashes -gt 0) { [void]$builder.Append([char]92, (2 * $slashes)) }
    [void]$builder.Append([char]34)
    return $builder.ToString()
}

function Get-AbsoluteNewPath {
    param([string]$Value)
    if ([IO.Path]::IsPathRooted($Value)) { return [IO.Path]::GetFullPath($Value) }
    return [IO.Path]::GetFullPath((Join-Path (Get-Location).Path $Value))
}

if (-not (Test-Path -LiteralPath $InputFile -PathType Leaf)) {
    throw "Input file not found: $InputFile"
}
if (-not (Test-Path -LiteralPath $Executable -PathType Leaf)) {
    throw "Executable not found: $Executable. Build the program first."
}
$inputPath = (Resolve-Path -LiteralPath $InputFile).ProviderPath
$exePath = (Resolve-Path -LiteralPath $Executable).ProviderPath
$outPath = Get-AbsoluteNewPath $OutputDirectory
$inputExtension = [IO.Path]::GetExtension($inputPath).ToLowerInvariant()
if ($inputExtension -eq '.gz') {
    throw 'Input must be uncompressed FASTQ/FASTA. Run Prepare-Inputs.ps1 first.'
}
if ($inputExtension -notin @('.fastq', '.fq', '.fasta', '.fa')) {
    throw 'Input extension must be .fastq, .fq, .fasta or .fa.'
}
if (Test-Path -LiteralPath $outPath) {
    throw "Output directory already exists: $outPath. Use a new directory; LightAlign appends to PAF files."
}
if ($UnitLength -lt $Window) { throw 'UnitLength must be at least Window.' }
[void][IO.Directory]::CreateDirectory($outPath)

$metricsPath = Join-Path $outPath 'kernel_metrics.json'
$stdoutPath = Join-Path $outPath 'stdout.log'
$stderrPath = Join-Path $outPath 'stderr.log'
$samplesPath = Join-Path $outPath 'memory_samples.csv'
$summaryJsonPath = Join-Path $outPath 'resource_summary.json'
$summaryCsvPath = Join-Path $outPath 'resource_summary.csv'
$pafPath = Join-Path $outPath 'overlap.paf'
$arguments = @('-i', $inputPath, '-O', $outPath,
    '-g', $GroupSize.ToString($invariant), '-w', $Window.ToString($invariant),
    '-e', $ErrorRate.ToString('R', $invariant), '-d', $DBA.ToString($invariant),
    '-l', $UnitLength.ToString($invariant), '-t', $Threads.ToString($invariant))
$argumentString = (($arguments | ForEach-Object { ConvertTo-NativeArgument ([string]$_) }) -join ' ')
[IO.File]::WriteAllText((Join-Path $outPath 'command.txt'),
    (ConvertTo-NativeArgument $exePath) + ' ' + $argumentString + [Environment]::NewLine, $utf8)

$info = New-Object Diagnostics.ProcessStartInfo
$info.FileName = $exePath
$info.Arguments = $argumentString
$info.WorkingDirectory = Split-Path -Parent $exePath
$info.UseShellExecute = $false
$info.CreateNoWindow = $true
$info.RedirectStandardOutput = $true
$info.RedirectStandardError = $true
$info.EnvironmentVariables['LIGHTALIGN_METRICS_FILE'] = $metricsPath
$process = New-Object Diagnostics.Process
$process.StartInfo = $info
$sampleWriter = $null
$stdoutStream = $null
$stderrStream = $null
$stdoutTask = $null
$stderrTask = $null
$watch = [Diagnostics.Stopwatch]::StartNew()
$startedAt = (Get-Date).ToString('o')
$processStarted = $false
$processId = $null
$exitCode = $null
$sampledPeak = [long]0
$sampleCount = 0
$status = 'not_started'
$failureMessage = $null
$finalMetrics = $null
$completionMarkerFound = $false
$reportedWorkerThreads = $null

try {
    $stdoutStream = New-Object IO.FileStream($stdoutPath, [IO.FileMode]::CreateNew,
        [IO.FileAccess]::Write, [IO.FileShare]::Read)
    $stderrStream = New-Object IO.FileStream($stderrPath, [IO.FileMode]::CreateNew,
        [IO.FileAccess]::Write, [IO.FileShare]::Read)
    $sampleWriter = New-Object IO.StreamWriter($samplesPath, $false, $utf8)
    $sampleWriter.AutoFlush = $true
    $sampleWriter.WriteLine('elapsed_seconds,working_set_bytes,peak_working_set_so_far_bytes,private_bytes,os_thread_count')
    if (-not $process.Start()) { throw 'Windows could not start LightAlign.' }
    $processStarted = $true
    $processId = $process.Id
    $status = 'running'
    [IO.File]::WriteAllText((Join-Path $outPath 'process_id.txt'), [string]$processId, $utf8)
    # Stream directly to files. No PowerShell callbacks and no full-log buffering.
    $stdoutTask = $process.StandardOutput.BaseStream.CopyToAsync($stdoutStream)
    $stderrTask = $process.StandardError.BaseStream.CopyToAsync($stderrStream)
    Write-Host "LightAlign started. PID: $processId; requested workers: $Threads"
    Write-Host "Output and live memory samples: $outPath"

    do {
        try {
            $process.Refresh()
            if (-not $process.HasExited) {
                $workingSet = [long]$process.WorkingSet64
                $peakSoFar = [long]$process.PeakWorkingSet64
                $privateBytes = [long]$process.PrivateMemorySize64
                $threadCount = $process.Threads.Count
                if ($peakSoFar -gt $sampledPeak) { $sampledPeak = $peakSoFar }
                $sampleWriter.WriteLine([string]::Join(',', @(
                    $watch.Elapsed.TotalSeconds.ToString('F6', $invariant),
                    $workingSet.ToString($invariant), $peakSoFar.ToString($invariant),
                    $privateBytes.ToString($invariant), $threadCount.ToString($invariant))))
                $sampleCount++
            }
        } catch [InvalidOperationException] {
            # A sample may race with process termination; final kernel metrics
            # are written by the executable and do not depend on these samples.
        }
        $finished = $process.WaitForExit($SampleIntervalMilliseconds)
    } while (-not $finished)
    $process.WaitForExit()
    $exitCode = $process.ExitCode
    [void]$stdoutTask.GetAwaiter().GetResult()
    [void]$stderrTask.GetAwaiter().GetResult()
    $stdoutStream.Dispose()
    $stderrStream.Dispose()
    $stdoutStream = $null
    $stderrStream = $null
    foreach ($line in [IO.File]::ReadLines($stdoutPath)) {
        if ($line.Contains('All steps done.')) { $completionMarkerFound = $true }
        if ($line -match 'Number of available threads:\s*(\d+)') {
            $reportedWorkerThreads = [int]$Matches[1]
        }
    }
    if ($exitCode -eq 0 -and (Test-Path -LiteralPath $pafPath -PathType Leaf) -and
        (Get-Item -LiteralPath $pafPath).Length -gt 0 -and $completionMarkerFound -and
        $reportedWorkerThreads -eq $Threads) {
        $status = 'succeeded'
    } else {
        $status = 'failed'
        $failureMessage = "Exit code: $exitCode; completion marker: $completionMarkerFound; reported workers: $reportedWorkerThreads. A nonempty overlap.paf and the requested thread count are required."
    }
} catch {
    $failureMessage = $_.Exception.Message
    if ($_.Exception -is [Management.Automation.PipelineStoppedException]) {
        $status = 'interrupted'
    } else {
        $status = 'runner_error'
    }
} finally {
    if ($processStarted) {
        try {
            if (-not $process.HasExited) {
                if ($status -eq 'running') { $status = 'interrupted' }
                $process.Kill()
                [void]$process.WaitForExit(10000)
            }
            if ($process.HasExited) { $exitCode = $process.ExitCode }
        } catch { }
    }
    foreach ($task in @($stdoutTask, $stderrTask)) {
        if ($null -ne $task) {
            try { [void]$task.Wait(10000) } catch { }
        }
    }
    foreach ($stream in @($sampleWriter, $stdoutStream, $stderrStream)) {
        if ($null -ne $stream) { try { $stream.Dispose() } catch { } }
    }
    $watch.Stop()
    if (Test-Path -LiteralPath $metricsPath -PathType Leaf) {
        try {
            $candidate = Get-Content -LiteralPath $metricsPath -Raw -Encoding UTF8 | ConvertFrom-Json
            if ($candidate.process_id -eq $processId) { $finalMetrics = $candidate }
        } catch { }
    }
    $authoritativePeak = $null
    $peakMiB = $null
    $peakCommit = $null
    $cpuUser = $null
    $cpuKernel = $null
    $inProcessElapsed = $null
    $metricsAvailable = $false
    if ($null -ne $finalMetrics) {
        $inProcessElapsed = $finalMetrics.elapsed_seconds
        $cpuUser = $finalMetrics.cpu_user_seconds
        $cpuKernel = $finalMetrics.cpu_kernel_seconds
        if ($finalMetrics.memory_query_succeeded) {
            $metricsAvailable = $true
            $authoritativePeak = [long]$finalMetrics.peak_working_set_bytes
            $peakMiB = [Math]::Round($authoritativePeak / 1MB, 3)
            $peakCommit = [long]$finalMetrics.peak_commit_bytes
        }
    }
    $summary = [PSCustomObject][ordered]@{
        status = $status
        exit_code = $exitCode
        process_id = $processId
        input_file = $inputPath
        output_directory = $outPath
        executable = $exePath
        requested_worker_threads = $Threads
        reported_worker_threads = $reportedWorkerThreads
        completion_marker_found = $completionMarkerFound
        group_size = $GroupSize
        window = $Window
        error_rate = $ErrorRate
        dba = $DBA
        alignment_unit_length = $UnitLength
        started_at_local = $startedAt
        finished_at_local = (Get-Date).ToString('o')
        launcher_elapsed_seconds = [Math]::Round($watch.Elapsed.TotalSeconds, 6)
        lightalign_elapsed_seconds = $inProcessElapsed
        cpu_user_seconds = $cpuUser
        cpu_kernel_seconds = $cpuKernel
        kernel_peak_available = $metricsAvailable
        peak_working_set_bytes = $authoritativePeak
        peak_working_set_mib = $peakMiB
        peak_commit_bytes = $peakCommit
        sampled_peak_working_set_bytes = $sampledPeak
        sample_count = $sampleCount
        memory_metric = 'Windows process peak working set; final in-process GetProcessMemoryInfo'
        paf_file = $pafPath
        error_message = $failureMessage
    }
    [IO.File]::WriteAllText($summaryJsonPath, ($summary | ConvertTo-Json -Depth 4), $utf8)
    $summary | Export-Csv -LiteralPath $summaryCsvPath -NoTypeInformation -Encoding UTF8
    if ($processStarted) { try { $process.Dispose() } catch { } }
}

Write-Host "Status: $status"
if ($null -ne $peakMiB) {
    Write-Host "Peak physical working set: $peakMiB MiB"
} else {
    Write-Warning 'Final kernel peak is unavailable; do not report the sampled peak as a completed-run peak.'
}
Write-Host "Resource records: $summaryCsvPath"
Write-Host "Program output: $stdoutPath"
Write-Host "Error output: $stderrPath"
if ($status -ne 'succeeded') {
    throw "Run did not succeed. See resource_summary.json and stderr.log in $outPath"
}
