[CmdletBinding(DefaultParameterSetName='Run')]
param(
    [Parameter(Mandatory,ParameterSetName='Run')][ValidateSet('CH','pooled')][string]$Case,
    [Parameter(Mandatory,ParameterSetName='Run')][string]$OutDir,
    [ValidatePattern('^[A-Za-z]\w*$')][string]$FunctionName = 'LL_mxl',
    [string]$BaselineDir = '',
    [string]$ReplicationRoot = 'C:\Users\miq\Documents\lasy\replication_package',
    [string]$Matlab = 'C:\Program Files\MATLAB\R2026b\bin\matlab.exe',
    [ValidateRange(3,1000)][int]$Repeats = 5,
    [ValidateRange(5,300)][int]$MinSeconds = 20,
    [ValidateRange(100,5000)][int]$SampleMilliseconds = 1000,
    [Parameter(Mandatory,ParameterSetName='SelfTest')][switch]$SelfTest
)
$ErrorActionPreference = 'Stop'

# Native counters avoid localized perf-counter names and identify processes by PID.
Add-Type -TypeDefinition @'
using System;
using System.ComponentModel;
using System.Runtime.InteropServices;
public static class MxlMemoryCounters {
    [StructLayout(LayoutKind.Sequential)]
    public struct Counters {
        public uint Size, PageFaults;
        public UIntPtr PeakWorkingSet, WorkingSet, QuotaPeakPagedPool, QuotaPagedPool;
        public UIntPtr QuotaPeakNonPagedPool, QuotaNonPagedPool, Pagefile, PeakPagefile, PrivateBytes;
    }
    [DllImport("kernel32.dll", SetLastError=true)]
    static extern IntPtr OpenProcess(uint access, bool inherit, int pid);
    [DllImport("kernel32.dll")] static extern bool CloseHandle(IntPtr handle);
    [DllImport("psapi.dll", SetLastError=true)]
    static extern bool GetProcessMemoryInfo(IntPtr handle, ref Counters counters, uint size);
    public static Counters Read(int pid) {
        IntPtr handle = OpenProcess(0x1000, false, pid);
        if (handle == IntPtr.Zero) throw new Win32Exception(Marshal.GetLastWin32Error());
        try {
            Counters c = new Counters();
            c.Size = (uint)Marshal.SizeOf(typeof(Counters));
            if (!GetProcessMemoryInfo(handle, ref c, c.Size))
                throw new Win32Exception(Marshal.GetLastWin32Error());
            return c;
        } finally { CloseHandle(handle); }
    }
}
'@

function Get-FaultDelta([double]$New,[double]$Old) {
    if ($New -ge $Old) { return $New-$Old }
    return $New-$Old+4294967296 # DWORD counter wrap.
}

if ($SelfTest) {
    $before = [MxlMemoryCounters]::Read($PID)
    $buffer = [byte[]]::new(1048576)
    for ($i=0;$i -lt $buffer.Length;$i+=4096) { $buffer[$i] = 1 }
    $after = [MxlMemoryCounters]::Read($PID)
    if ($after.WorkingSet.ToUInt64() -eq 0 -or $after.PrivateBytes.ToUInt64() -eq 0 -or
        (Get-FaultDelta $after.PageFaults $before.PageFaults) -lt 0 -or
        (Get-FaultDelta 3 4294967294) -ne 5) { throw 'Counter self-check failed.' }
    Write-Output 'Native process memory/page-fault counter self-check passed.'
    return
}

$outPath = [IO.Path]::GetFullPath($OutDir).TrimEnd('\')
$frozenPath = [IO.Path]::GetFullPath($ReplicationRoot).TrimEnd('\')
if ($outPath -eq $frozenPath -or $outPath.StartsWith($frozenPath+'\',[StringComparison]::OrdinalIgnoreCase)) {
    throw 'Benchmark outputs must be outside replication_package.'
}
if (Test-Path -LiteralPath $outPath) { throw 'Use a new output directory for each benchmark series.' }
if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) {
    throw 'A MATLAB session already exists. Confirm the computer is free before launching this batch.'
}
New-Item -ItemType Directory -Path $outPath | Out-Null
function Quote-Matlab([string]$Text) { return "'"+$Text.Replace("'","''")+"'" }
$quotedArgs = (@($Case,$FunctionName,$outPath,$ReplicationRoot,$BaselineDir) | ForEach-Object { Quote-Matlab $_ }) -join ','
$expression = 'addpath('+ (Quote-Matlab $PSScriptRoot) +'); bench_mxl_memory('+$quotedArgs+','+$Repeats+','+$MinSeconds+');'
$arguments = '-wait -batch "'+$expression+'" -logfile "'+(Join-Path $outPath 'matlab.log')+'"'
$matlabProcess = Start-Process -FilePath $Matlab -ArgumentList $arguments -WindowStyle Hidden -PassThru
$readyPath = Join-Path $outPath 'ready.json'
$resultPath = Join-Path $outPath 'evaluation.json'
$rows = [Collections.Generic.List[object]]::new()
$previous = @{}
$timer = [Diagnostics.Stopwatch]::StartNew()

function Sample-Processes {
    $sampleNumber = $rows.Count / $processIds.Count
    foreach ($processId in $processIds) {
        $clock = $timer.Elapsed.TotalSeconds
        $c = [MxlMemoryCounters]::Read($processId)
        $delta = $null; $rate = $null
        if ($previous.ContainsKey($processId)) {
            $delta = Get-FaultDelta $c.PageFaults $previous[$processId].Faults
            $rate = $delta / ($clock-$previous[$processId].Seconds)
        }
        $rows.Add([pscustomobject]@{
            TimestampUTC = [DateTime]::UtcNow.ToString('o'); ElapsedSeconds = $clock; Sample = $sampleNumber
            PID = $processId; Role = $(if ($processId -eq $processIds[0]) { 'client' } else { 'worker'+[array]::IndexOf($processIds,$processId) })
            PageFaultCount = $c.PageFaults; PageFaultDelta = $delta; PageFaultsPerSecond = $rate
            WorkingSetBytes = $c.WorkingSet.ToUInt64(); PrivateBytes = $c.PrivateBytes.ToUInt64()
            LifetimePeakWorkingSetBytes = $c.PeakWorkingSet.ToUInt64()
            LifetimePeakPrivateCommitBytes = $c.PeakPagefile.ToUInt64()
        })
        $previous[$processId] = @{ Faults=$c.PageFaults; Seconds=$clock }
    }
}

try {
    while (-not (Test-Path -LiteralPath $readyPath)) {
        $matlabProcess.Refresh()
        if ($matlabProcess.HasExited) { throw "MATLAB exited before measurement (code $($matlabProcess.ExitCode)); inspect matlab.log." }
        if ($timer.Elapsed.TotalSeconds -gt 600) { throw 'Pool/data/warm-up exceeded the sampler startup timeout.' }
        Start-Sleep -Milliseconds 250
    }
    $ready = Get-Content -LiteralPath $readyPath -Raw | ConvertFrom-Json
    $processIds = @($ready.pids | ForEach-Object { [int]$_ })
    if ($processIds.Count -ne 4 -or @($processIds | Select-Object -Unique).Count -ne 4) {
        throw 'Expected one distinct client and three distinct worker PIDs.'
    }
    Sample-Processes
    New-Item -ItemType File -Path (Join-Path $outPath 'sample.start') | Out-Null
    while (-not (Test-Path -LiteralPath $resultPath)) {
        Start-Sleep -Milliseconds $SampleMilliseconds
        $matlabProcess.Refresh()
        if ($matlabProcess.HasExited) { throw "MATLAB exited during measurement (code $($matlabProcess.ExitCode)); inspect matlab.log." }
        Sample-Processes
    }
    $rows | Export-Csv -LiteralPath (Join-Path $outPath 'process_samples.tsv') -Delimiter "`t" -NoTypeInformation
    $perProcess = foreach ($processId in $processIds) {
        $samples = @($rows | Where-Object PID -eq $processId)
        $first = $samples[0]; $last = $samples[-1]
        [pscustomobject]@{
            PID = $processId; Role = $first.Role; SampleSeconds = $last.ElapsedSeconds-$first.ElapsedSeconds
            PageFaults = ($samples | Measure-Object -Property PageFaultDelta -Sum).Sum
            MeanPageFaultsPerSecond = ($samples | Measure-Object -Property PageFaultDelta -Sum).Sum / ($last.ElapsedSeconds-$first.ElapsedSeconds)
            PeakPageFaultsPerSecond = ($samples | Measure-Object -Property PageFaultsPerSecond -Maximum).Maximum
            SampledPeakWorkingSetMiB = ($samples | Measure-Object -Property WorkingSetBytes -Maximum).Maximum / 1MB
            SampledPeakPrivateMiB = ($samples | Measure-Object -Property PrivateBytes -Maximum).Maximum / 1MB
            LifetimePeakWorkingSetMiB = $last.LifetimePeakWorkingSetBytes / 1MB
            LifetimePeakPrivateCommitMiB = $last.LifetimePeakPrivateCommitBytes / 1MB
        }
    }
    $totals = @($rows | Group-Object Sample | ForEach-Object {
        [pscustomobject]@{
            WorkingSet = ($_.Group | Measure-Object -Property WorkingSetBytes -Sum).Sum
            Private = ($_.Group | Measure-Object -Property PrivateBytes -Sum).Sum
        }
    })
    $result = Get-Content -LiteralPath $resultPath -Raw | ConvertFrom-Json
    [pscustomobject]@{
        Evaluation = $result; ProcessMetrics = @($perProcess); SampleMilliseconds = $SampleMilliseconds
        SampledPeakTotalWorkingSetMiB = ($totals | Measure-Object -Property WorkingSet -Maximum).Maximum / 1MB
        SampledPeakTotalPrivateMiB = ($totals | Measure-Object -Property Private -Maximum).Maximum / 1MB
        WarmupEvaluations = 2; CounterSource = 'GetProcessMemoryInfo; cumulative DWORD PageFaultCount'
        MemoryPeakScope = 'Sampled peaks exclude warm-up; native lifetime peaks include pool startup and warm-up.'
    } | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $outPath 'summary.json') -Encoding UTF8
    $perProcess | Export-Csv -LiteralPath (Join-Path $outPath 'process_summary.tsv') -Delimiter "`t" -NoTypeInformation
    New-Item -ItemType File -Path (Join-Path $outPath 'sample.done') | Out-Null
    $matlabProcess.WaitForExit()
    if ($matlabProcess.ExitCode -ne 0) { throw "MATLAB exited with code $($matlabProcess.ExitCode); inspect matlab.log." }
    Write-Output (Join-Path $outPath 'summary.json')
} finally {
    if ($rows.Count -gt 0) {
        $rows | Export-Csv -LiteralPath (Join-Path $outPath 'process_samples.tsv') -Delimiter "`t" -NoTypeInformation
    }
    # Let the MATLAB handshake timeout close its own pool after sampler failure.
    $matlabProcess.Refresh()
    if (-not $matlabProcess.HasExited) { $matlabProcess.WaitForExit() }
}
