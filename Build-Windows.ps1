[CmdletBinding()]
param(
    [string]$CompilerPath
)

# Windows PowerShell 5.1 compatible. This script does not install software.
Set-StrictMode -Version 2.0
$ErrorActionPreference = 'Stop'
$package = $PSScriptRoot
$src = Join-Path $package 'src'
$bin = Join-Path $package 'bin'
$buildRoot = Join-Path $package 'build'
$stamp = (Get-Date).ToString('yyyyMMdd_HHmmss') + '_' + [Guid]::NewGuid().ToString('N').Substring(0, 8)
$work = Join-Path $buildRoot $stamp
[void][IO.Directory]::CreateDirectory($work)
[void][IO.Directory]::CreateDirectory($bin)
$utf8 = New-Object Text.UTF8Encoding($false)
$logPath = Join-Path $work 'build.log'
$started = (Get-Date).ToString('o')
[IO.File]::WriteAllText($logPath, "Build started: $started`r`nPackage: $package`r`n", $utf8)

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

function Invoke-NativeCaptured {
    param([string]$File, [string[]]$Arguments, [string]$WorkingDirectory)
    $info = New-Object Diagnostics.ProcessStartInfo
    $info.FileName = $File
    $info.Arguments = (($Arguments | ForEach-Object { ConvertTo-NativeArgument ([string]$_) }) -join ' ')
    $info.WorkingDirectory = $WorkingDirectory
    $info.UseShellExecute = $false
    $info.CreateNoWindow = $true
    $info.RedirectStandardOutput = $true
    $info.RedirectStandardError = $true
    $process = New-Object Diagnostics.Process
    $process.StartInfo = $info
    try {
        if (-not $process.Start()) { throw "Could not start $File" }
        $stdoutTask = $process.StandardOutput.ReadToEndAsync()
        $stderrTask = $process.StandardError.ReadToEndAsync()
        $process.WaitForExit()
        $stdout = $stdoutTask.GetAwaiter().GetResult()
        $stderr = $stderrTask.GetAwaiter().GetResult()
        return [PSCustomObject]@{ ExitCode = $process.ExitCode; Stdout = $stdout; Stderr = $stderr }
    } finally {
        $process.Dispose()
    }
}

function Find-MSVC {
    $cl = Get-Command cl.exe -ErrorAction SilentlyContinue
    if ($null -ne $cl -and ($env:VSCMD_ARG_TGT_ARCH -eq 'x64' -or [string]::IsNullOrEmpty($env:VSCMD_ARG_TGT_ARCH))) {
        return $cl.Source
    }
    $vswhere = $null
    foreach ($programFiles in @(${env:ProgramFiles(x86)}, $env:ProgramFiles)) {
        if (-not [string]::IsNullOrEmpty($programFiles)) {
            $candidate = Join-Path $programFiles 'Microsoft Visual Studio\Installer\vswhere.exe'
            if (Test-Path -LiteralPath $candidate -PathType Leaf) { $vswhere = $candidate; break }
        }
    }
    if ($null -eq $vswhere) {
        throw 'MSVC C++ compiler not found. Install Visual Studio Build Tools with Desktop development with C++, or specify -CompilerPath to a portable zig.exe.'
    }
    $query = Invoke-NativeCaptured -File $vswhere -Arguments @('-latest', '-products', '*', '-requires', 'Microsoft.VisualStudio.Component.VC.Tools.x86.x64', '-property', 'installationPath') -WorkingDirectory $work
    $installation = $query.Stdout.Trim()
    if ($query.ExitCode -ne 0 -or [string]::IsNullOrWhiteSpace($installation)) {
        throw 'Visual Studio C++ x64 build tools are not installed. Install the Desktop development with C++ workload, including MSVC and Windows SDK, then retry.'
    }
    $devcmd = Join-Path $installation 'Common7\Tools\VsDevCmd.bat'
    if (-not (Test-Path -LiteralPath $devcmd -PathType Leaf)) { throw "Missing developer environment script: $devcmd" }
    $info = New-Object Diagnostics.ProcessStartInfo
    $info.FileName = $env:ComSpec
    # cmd.exe requires its own outer quoting around a quoted batch filename.
    $info.Arguments = '/d /s /c ""' + $devcmd + '" -no_logo -arch=x64 -host_arch=x64 >nul && set"'
    $info.UseShellExecute = $false
    $info.CreateNoWindow = $true
    $info.RedirectStandardOutput = $true
    $info.RedirectStandardError = $true
    $process = New-Object Diagnostics.Process
    $process.StartInfo = $info
    try {
        [void]$process.Start()
        $outTask = $process.StandardOutput.ReadToEndAsync()
        $errTask = $process.StandardError.ReadToEndAsync()
        $process.WaitForExit()
        $environmentText = $outTask.GetAwaiter().GetResult()
        $environmentError = $errTask.GetAwaiter().GetResult()
        if ($process.ExitCode -ne 0) { throw "Could not initialize MSVC: $environmentError" }
        foreach ($line in ($environmentText -split "`r?`n")) {
            $equals = $line.IndexOf('=')
            if ($equals -gt 0) {
                [Environment]::SetEnvironmentVariable($line.Substring(0, $equals), $line.Substring($equals + 1), 'Process')
            }
        }
    } finally { $process.Dispose() }
    $cl = Get-Command cl.exe -ErrorAction SilentlyContinue
    if ($null -eq $cl) { throw 'MSVC environment initialized but cl.exe was not found.' }
    return $cl.Source
}

$oldGlobalCache = $env:ZIG_GLOBAL_CACHE_DIR
$oldLocalCache = $env:ZIG_LOCAL_CACHE_DIR
$oldTemp = $env:TEMP
$oldTmp = $env:TMP
try {
    $temporary = Join-Path $work 'tmp'
    [void][IO.Directory]::CreateDirectory($temporary)
    $env:TEMP = $temporary
    $env:TMP = $temporary
    $sourceNames = @('help.cpp', 'step4_5.cpp', 'step3.cpp', 'step2.cpp', 'step1.cpp', 'functions.cpp', 'MurmurHash3.cpp', 'preprocess_fasta.cpp')
    $sources = @()
    foreach ($name in $sourceNames) {
        $path = Join-Path $src $name
        if (-not (Test-Path -LiteralPath $path -PathType Leaf)) { throw "Source file missing: $path" }
        $sources += $path
    }
    if ([string]::IsNullOrWhiteSpace($CompilerPath)) {
        $compiler = Find-MSVC
    } else {
        if (-not (Test-Path -LiteralPath $CompilerPath -PathType Leaf)) { throw "Compiler not found: $CompilerPath" }
        $compiler = (Resolve-Path -LiteralPath $CompilerPath).ProviderPath
    }
    $compilerName = [IO.Path]::GetFileName($compiler).ToLowerInvariant()
    $stagedExe = Join-Path $work 'LightAlign.exe'
    if ($compilerName -eq 'zig.exe' -or $compilerName -eq 'zig') {
        $kind = 'Zig C++'
        $env:ZIG_GLOBAL_CACHE_DIR = Join-Path $buildRoot 'zig-cache-global'
        $env:ZIG_LOCAL_CACHE_DIR = Join-Path $buildRoot 'zig-cache-local'
        [void][IO.Directory]::CreateDirectory($env:ZIG_GLOBAL_CACHE_DIR)
        [void][IO.Directory]::CreateDirectory($env:ZIG_LOCAL_CACHE_DIR)
        $version = Invoke-NativeCaptured -File $compiler -Arguments @('version') -WorkingDirectory $work
        $arguments = @('c++', '-O2', '-std=c++17', '-fno-lto', '-D_FILE_OFFSET_BITS=64', '-DNOMINMAX', '-DWIN32_LEAN_AND_MEAN') + $sources + @('-o', $stagedExe, '-lpsapi')
    } elseif ($compilerName -eq 'cl.exe' -or $compilerName -eq 'cl') {
        $kind = 'MSVC C++'
        $version = Invoke-NativeCaptured -File $compiler -Arguments @('/Bv') -WorkingDirectory $work
        $arguments = @('/nologo', '/O2', '/std:c++17', '/EHsc', '/utf-8', '/MT', '/DNOMINMAX', '/DWIN32_LEAN_AND_MEAN', '/D_WIN32_WINNT=0x0601', ('/Fe' + $stagedExe), ('/Fo' + $work + '\')) + $sources + @('/link', '/MACHINE:X64', '/INCREMENTAL:NO', 'psapi.lib')
    } else {
        throw 'CompilerPath must point to cl.exe or zig.exe.'
    }
    $versionText = ($version.Stdout + $version.Stderr).Trim()
    [IO.File]::WriteAllText((Join-Path $work 'compiler_version.txt'), $versionText + [Environment]::NewLine, $utf8)
    $command = (ConvertTo-NativeArgument $compiler) + ' ' + (($arguments | ForEach-Object { ConvertTo-NativeArgument ([string]$_) }) -join ' ')
    [IO.File]::AppendAllText($logPath, "Compiler: $compiler`r`n$versionText`r`nCommand: $command`r`n", $utf8)
    Write-Host "Compiling with $kind ..."
    Write-Host "Build directory: $work"
    $result = Invoke-NativeCaptured -File $compiler -Arguments $arguments -WorkingDirectory $work
    [IO.File]::WriteAllText((Join-Path $work 'build.stdout.log'), $result.Stdout, $utf8)
    [IO.File]::WriteAllText((Join-Path $work 'build.stderr.log'), $result.Stderr, $utf8)
    [IO.File]::AppendAllText($logPath, "`r`nSTDOUT:`r`n" + $result.Stdout + "`r`nSTDERR:`r`n" + $result.Stderr + "`r`nExit code: " + $result.ExitCode + "`r`n", $utf8)
    if ($result.ExitCode -ne 0 -or -not (Test-Path -LiteralPath $stagedExe -PathType Leaf)) {
        throw "Compilation failed (exit $($result.ExitCode)). See $logPath"
    }
    $helpResult = Invoke-NativeCaptured -File $stagedExe -Arguments @('-h') -WorkingDirectory $work
    [IO.File]::WriteAllText((Join-Path $work 'help.txt'), $helpResult.Stdout + $helpResult.Stderr, $utf8)
    if ($helpResult.ExitCode -ne 0 -or $helpResult.Stdout -notmatch '\-t\s') {
        throw "Compiled program did not pass the help/threads check. See $work\help.txt"
    }
    $target = Join-Path $bin 'LightAlign.exe'
    if (Test-Path -LiteralPath $target -PathType Leaf) {
        Copy-Item -LiteralPath $target -Destination (Join-Path $work 'previous_LightAlign.exe')
    }
    Copy-Item -LiteralPath $stagedExe -Destination $target -Force
    $hashes = @(Get-ChildItem -LiteralPath $src -File | Where-Object { $_.Extension -in @('.cpp', '.h') } | ForEach-Object {
        [PSCustomObject]@{ file = $_.Name; sha256 = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLowerInvariant() }
    })
    $hashes | Export-Csv -LiteralPath (Join-Path $work 'source_hashes.csv') -NoTypeInformation -Encoding UTF8
    $manifest = [PSCustomObject][ordered]@{
        started_at_local = $started
        finished_at_local = (Get-Date).ToString('o')
        compiler_kind = $kind
        compiler_path = $compiler
        compiler_version = $versionText
        command = $command
        build_directory = $work
        executable = $target
        executable_sha256 = (Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash.ToLowerInvariant()
        source_files = $hashes
        help_check_passed = $true
    }
    [IO.File]::WriteAllText((Join-Path $work 'build_manifest.json'), ($manifest | ConvertTo-Json -Depth 6), $utf8)
    Write-Host "Build succeeded: $target"
    Write-Host "Build log and provenance: $work"
} catch {
    [IO.File]::AppendAllText($logPath, "`r`nBUILD FAILED: " + $_.Exception.Message + "`r`n", $utf8)
    Write-Host "Build log: $logPath"
    throw
} finally {
    [Environment]::SetEnvironmentVariable('ZIG_GLOBAL_CACHE_DIR', $oldGlobalCache, 'Process')
    [Environment]::SetEnvironmentVariable('ZIG_LOCAL_CACHE_DIR', $oldLocalCache, 'Process')
    [Environment]::SetEnvironmentVariable('TEMP', $oldTemp, 'Process')
    [Environment]::SetEnvironmentVariable('TMP', $oldTmp, 'Process')
}
