[CmdletBinding(PositionalBinding = $false)]
param(
    [string]$Python,
    [string]$ComfyRoot,
    [string]$FrontendRoot,
    [Alias("h")][switch]$Help,
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$ComfyArgs
)

$ErrorActionPreference = "Stop"

if ($Help) {
    Write-Output @"
Usage: powershell -File scripts/start.ps1 [-Python <python.exe>] [-ComfyRoot <path>] [-FrontendRoot <path>] [ComfyUI arguments...]

Starts ComfyUI with the sibling OpenBio frontend dist, local assets,
and API nodes disabled. Additional arguments are forwarded to ComfyUI.

ComfyUI can also be selected with OPENBIO_COMFYUI_ROOT. The frontend repository
can also be selected with OPENBIO_FRONTEND_ROOT.
"@
    exit 0
}

$PluginRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot "resolve_runtime.ps1")

function Resolve-OpenBioFrontendRoot([string]$Requested) {
    if ($Requested) {
        $Candidates = @($Requested)
        $Source = "-FrontendRoot"
    } elseif ($env:OPENBIO_FRONTEND_ROOT) {
        $Candidates = @($env:OPENBIO_FRONTEND_ROOT)
        $Source = "OPENBIO_FRONTEND_ROOT"
    } else {
        $Candidates = @(
            (Join-Path (Split-Path -Parent $PluginRoot) "ComfyUI_frontend"),
            (Join-Path (Split-Path -Parent $ComfyRoot) "ComfyUI_frontend")
        )
        $Source = $null
    }

    foreach ($Candidate in $Candidates) {
        if (Test-Path -LiteralPath $Candidate -PathType Container) {
            return (Resolve-Path -LiteralPath $Candidate).Path
        }
    }

    if ($Source) {
        throw "$Source does not point to a frontend repository: $($Candidates[0])"
    }
    throw "ComfyUI frontend repository was not found. Pass -FrontendRoot or set OPENBIO_FRONTEND_ROOT."
}

$ComfyRoot = Resolve-OpenBioComfyRoot $ComfyRoot $PluginRoot
$FrontendRoot = Resolve-OpenBioFrontendRoot $FrontendRoot
$FrontendDist = Join-Path $FrontendRoot "dist"

$FrontendIndex = Join-Path $FrontendDist "index.html"
if ((-not (Test-Path -LiteralPath $FrontendIndex -PathType Leaf)) -or ((Get-Item -LiteralPath $FrontendIndex).Length -eq 0)) {
    throw @"
OpenBio frontend dist is incomplete: $FrontendDist
Missing or empty: index.html
Build it first:
  Set-Location "$FrontendRoot"
  corepack pnpm install --frozen-lockfile
  corepack pnpm build:openbio
"@
}

$PythonExe = Resolve-OpenBioPython $Python $ComfyRoot
& $PythonExe -c "import sys; raise SystemExit(sys.version_info < (3, 12))"
if ($LASTEXITCODE -ne 0) {
    throw "openbio-singlecell requires Python 3.12 or newer."
}

$Arguments = @(
    (Join-Path $ComfyRoot "main.py"),
    "--disable-api-nodes",
    "--enable-assets",
    "--front-end-root",
    $FrontendDist
)
$HasCacheMode = @($ComfyArgs | Where-Object { $_ -match '^--(cache-(classic|none|lru|ram)|high-ram)(=|$)' }).Count -gt 0
if (-not $HasCacheMode) {
    $Arguments += "--cache-classic"
}
$Arguments += $ComfyArgs

Write-Host "Using Python: $PythonExe"
Push-Location $ComfyRoot
try {
    & $PythonExe @Arguments
    $ExitCode = $LASTEXITCODE
} finally {
    Pop-Location
}
exit $ExitCode
